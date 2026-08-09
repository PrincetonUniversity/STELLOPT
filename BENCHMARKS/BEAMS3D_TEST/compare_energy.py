#!/usr/bin/env python3
# Energy-conservation check for BEAMS3D full-orbit (RHO_FULLORBIT) runs.
#
# With no electrostatic potential the Lorentz equation conserves |v| exactly for
# an arbitrary B(x), because the v x B force is perpendicular to v.  Any drift is
# therefore integration error or a bug, which makes this a strong regression test
# that is insensitive to the random gyro-phase used when stepping from the gyro
# centre to the particle position.
#
# The speed is rebuilt from the parallel velocity, the magnetic moment and |B|,
#
#     v^2 = vll^2 + 2 * moment * |B| / m
#
# rather than from the vr/vphi/vz that out_beams3d_part also writes.  The two
# agree to machine precision (4.6E-16 relative on the benchmark cases); this
# route is used because vll/moment/B_lines are present in every BEAMS3D output
# while the velocity vector is not.
#
# Unlike compare.py this script exits non-zero on failure: the threshold is
# deterministic, so a failure here is a real regression rather than statistics.
import sys, os
import numpy as np
from argparse import ArgumentParser
sys.path.insert(0, '../../pySTEL/')
from libstell.beams3d import BEAMS3D

_E_CHARGE_ = 1.602176634E-19

# Main routine
if __name__=="__main__":
    # Parse Arguments
    parser = ArgumentParser(description=
        '''Checks energy conservation along BEAMS3D full-orbit trajectories.''')
    parser.add_argument("--file", dest="run_name",
        help="File to process", default = None)
    parser.add_argument("--tol", dest="failtol", type=float,
        help="Largest tolerated |dE/E| along a trajectory", default = 1.0E-5)

    args = parser.parse_args()
    run_name = args.run_name

    # Do nothing if no filename
    if not run_name: sys.exit(0)
    failtol = args.failtol
    b3d = BEAMS3D()
    try:
        b3d.read_beams3d(f'beams3d_{run_name}.h5')
    except:
        print(f'  ERROR: Cannot find run {run_name}')
        sys.exit(-1)

    # Rebuild the speed from vll, mu and |B|.  Index 0 holds the launch point,
    # which is a gyro-centre quantity and is skipped.
    vll  = b3d.vll_lines[1:,:]
    mu   = b3d.moment_lines[1:,:]
    modb = b3d.B_lines[1:,:]
    rr   = b3d.R_lines[1:,:]
    mass = b3d.mass[None,:]
    energy = 0.5*mass*(vll*vll + 2.0*mu*modb/mass)/_E_CHARGE_
    valid  = (rr > 0) & (modb > 0) & np.isfinite(energy)

    print(f'BEAMS3D VERSION: {b3d.VERSION:4.2f}')
    print('=================')
    print(f'  Full orbit energy conservation, tol = {failtol:9.3e}')
    print(f'  Marker -- E0 [keV] -- max|dE/E| -- dE/E(end)')
    lfail = False
    for i in range(energy.shape[1]):
        trace = energy[valid[:,i],i]
        if trace.size < 2:
            # A marker with no usable points means the trajectory never ran or
            # went non-finite (e.g. INT_TYPE='RKH68', which silently produces
            # NaN in full orbit).  That must not be reported as a pass.
            print(f'  {i} NO USABLE POINTS')
            lfail = True
            continue
        e0    = trace[0]
        dmax  = np.max(np.abs(trace-e0))/e0
        dend  = (trace[-1]-e0)/e0
        print(f'  {i} {e0*1.0E-3:9.4f} {dmax:11.4e} {dend:+11.4e}')
        if dmax > failtol:
            lfail = True
    print('=================')

    # Error Status
    if lfail:
        print('  STATUS: FAIL!!!!!')
        sys.exit(1)
    else:
        print('  STATUS: PASS')
