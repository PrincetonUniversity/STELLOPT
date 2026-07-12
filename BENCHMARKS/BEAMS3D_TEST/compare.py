#!/usr/bin/env python3
import sys, os
import json
import h5py
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
        lost = b3d.end_state == 2
        valid = b3d.wall_hit_valid == 1
        assert b3d.lplasma_only
        assert b3d.lwall_from_vmec
        assert np.array_equal(valid, lost)
        assert np.array_equal(b3d.wall_hit_field_valid == 1, valid)
        assert np.all(b3d.wall_hit_model[valid] == 1)
        assert np.all(b3d.wall_hit_model[~valid] == 0)
        assert np.all((b3d.wall_hit_fraction[valid] >= 0.0) &
                      (b3d.wall_hit_fraction[valid] <= 1.0))
        assert np.all(b3d.wall_hit_time[valid] >= 0.0)
        assert np.all(b3d.wall_hit_time[valid] <= b3d.t_end[valid])
        assert np.all(np.isfinite(b3d.wall_hit_r[valid]))
        assert np.all(np.isfinite(b3d.wall_hit_phi[valid]))
        assert np.all(np.isfinite(b3d.wall_hit_z[valid]))
        for values in (b3d.wall_hit_vll, b3d.wall_hit_moment, b3d.wall_hit_b,
                       b3d.wall_hit_s, b3d.wall_hit_u, b3d.wall_hit_energy):
            assert np.all(np.isfinite(values[valid]))
            assert np.all(values[~valid] < -1.0e300)
        assert np.all(b3d.wall_hit_b[valid] > 0.0)
        assert np.all(b3d.wall_hit_energy[valid] > 0.0)
        expected_energy = (0.5*b3d.mass[valid]*b3d.wall_hit_vll[valid]**2 +
                           b3d.wall_hit_moment[valid]*b3d.wall_hit_b[valid])
        assert np.array_equal(b3d.wall_hit_energy[valid], expected_energy)
        assert np.all(b3d.wall_hit_face[~valid] == -1)
        assert np.all(b3d.wall_hit_fraction[~valid] == -1.0)
        assert np.all(b3d.wall_hit_time[~valid] == -1.0)
        assert np.all(b3d.wall_hit_r[~valid] < -1.0e300)
        assert np.all(b3d.wall_hit_phi[~valid] < -1.0e300)
        assert np.all(b3d.wall_hit_z[~valid] < -1.0e300)
        with h5py.File(f'beams3d_{run_name}.h5', 'r') as event_data:
            r_lines = event_data['R_lines'][:]
            z_lines = event_data['Z_lines'][:]
            phi_lines = event_data['PHI_lines'][:]
            vll_lines = event_data['vll_lines'][:]
            moment_lines = event_data['moment_lines'][:]
            b_lines = event_data['B_lines'][:]
            s_lines = event_data['S_lines'][:]
            u_lines = event_data['U_lines'][:]
            time_lines = event_data['time_lines'][:]
            wall_strikes = event_data['wall_strikes'][:]
            for name in ('wall_hit_valid', 'wall_hit_field_valid', 'wall_hit_model',
                         'wall_hit_face', 'wall_hit_fraction',
                         'wall_hit_time', 'wall_hit_r', 'wall_hit_phi',
                         'wall_hit_z', 'wall_hit_vll', 'wall_hit_moment',
                         'wall_hit_b', 'wall_hit_s', 'wall_hit_u',
                         'wall_hit_energy', 'time_lines'):
                assert 'description' in event_data[name].attrs
        populated = (r_lines[valid] != 0.0) | (z_lines[valid] != 0.0)
        all_populated = (r_lines != 0.0) | (z_lines != 0.0)
        assert np.all(time_lines[all_populated] >= 0.0)
        for particle in range(time_lines.shape[0]):
            particle_time = time_lines[particle, all_populated[particle]]
            assert np.all(np.diff(particle_time) >= 0.0)
        last = populated.shape[1] - 1 - np.argmax(populated[:, ::-1], axis=1)
        rows = np.flatnonzero(valid)
        assert np.sum(wall_strikes) == rows.size
        expected_strikes = np.bincount(b3d.wall_hit_face[valid] - 1,
                                       minlength=wall_strikes.size)
        assert np.array_equal(wall_strikes, expected_strikes)
        assert np.array_equal(b3d.wall_hit_r[valid], r_lines[rows, last])
        assert np.array_equal(b3d.wall_hit_phi[valid], phi_lines[rows, last])
        assert np.array_equal(b3d.wall_hit_z[valid], z_lines[rows, last])
        assert np.array_equal(b3d.wall_hit_vll[valid], vll_lines[rows, last])
        assert np.array_equal(b3d.wall_hit_moment[valid], moment_lines[rows, last])
        assert np.array_equal(b3d.wall_hit_b[valid], b_lines[rows, last])
        assert np.array_equal(b3d.wall_hit_s[valid], s_lines[rows, last])
        assert np.array_equal(b3d.wall_hit_u[valid], u_lines[rows, last])
        assert np.array_equal(b3d.wall_hit_time[valid], time_lines[rows, last])
    if run_name in ['ORBITS_slow']:
        data['vll0']=b3d.vll_lines[0,:].flatten().tolist()
        data['epower']=(b3d.epower_prof[0,:].flatten()*1.0E20).tolist() # factor to make values visible in output.
        data['ipower']=(b3d.ipower_prof[0,:].flatten()*1.0E20).tolist()
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
        sys.exit(1)
    else:
        print('  STATUS: PASS')
        sys.exit(0)
