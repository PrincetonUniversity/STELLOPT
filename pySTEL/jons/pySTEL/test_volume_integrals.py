# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 14:10:24 2017

@author: jons
"""

import sys
import numpy as np
import matplotlib.pyplot as plt


def import_path(p):
    # check if p is already in sys.path, append if not
    try:
        sys.path.index(p)
    except ValueError:
        sys.path.append(p)

#MHDpythonPath="Z:\\03_Masterarbeit\\00_programs\\MHDpython\\"
MHDpythonPath="/home/jonathan/Uni/03_Masterarbeit/00_programs/MHDpython/"
import_path(MHDpythonPath+"pyVMEC")


from read_vmec import read_vmec

extension="w7x_ref_84"
#datadir="Z:\\03_Masterarbeit\\01_analysis\\24_magnetic_coordinates\\"
datadir="/home/jonathan/Uni/03_Masterarbeit/data/w7x_ref_84/"

do_plot_dVoldPsi = False
do_plot_iota_profile = False


# construct vmec filename
vmec_filename="wout."+extension+".nc"

# read vmec data
vmec_data=read_vmec(datadir+vmec_filename)

# number of flux surfaces
ns=vmec_data['ns']

# total LCFS volume
vol=vmec_data['Volume']

# Volume derivative (half mesh)
vp=np.squeeze(vmec_data['vp'])

# toroidal flux (full mesh)
phi=np.squeeze(vmec_data['phi'])

# phipf = phi/(sfull*sfull)   (full mesh)
phipf=np.squeeze(vmec_data['phipf'])

# s grids: full and half
sfull=np.sqrt(np.linspace(0,ns-1,ns)/(ns-1))
shalf=np.sqrt(np.linspace(0.5,ns-0.5,ns)/(ns-1))

# four pi squared
fps = (2.0*np.pi)*(2.0*np.pi)

# as calculated by RefEq webservice
volume_profile=np.cumsum(vp*fps/(ns-1))  # => gives correct total LCFS volume in last point
print("integration precision: %.3e m^3" %np.abs(volume_profile[-1] - vol))

plt.figure()
plt.plot(1.0, vol, 'ko')   # target point
plt.plot(sfull*sfull, volume_profile, 'k')
plt.xlabel("normalized radius")
plt.ylabel("enclosed volume / m^3")
plt.title(extension)

#mean_vol_deriv = np.mean(volume_profile[1:]/(sfull[1:]*sfull[1:]))
#plt.figure()
#plt.plot(sfull*sfull, volume_profile-(sfull*sfull)*mean_vol_deriv, 'k')
#plt.xlabel("normalized radius")
#plt.ylabel("variation of enclosed volume / m^3")
#plt.title(extension)


#plt.figure()
#plt.plot(sfull*sfull, phi)
#plt.xlabel("normalized radius")
#plt.ylabel("toroidal flux / Wb")
#plt.title(extension)


## => numerical noise
#mean_flux_deriv = np.mean(phi[1:]/(sfull[1:]*sfull[1:]))
#plt.figure()
#plt.plot(sfull*sfull, phi-(sfull*sfull)*mean_flux_deriv)
#plt.xlabel("normalized radius")
#plt.ylabel("variation of toroidal flux / Wb")
#plt.title(extension)



if do_plot_dVoldPsi:
    # => so really all non-proportional variation here stems from the volume derivative
    plt.figure()
    # omit magnetic axis => no volume there
    plt.plot(sfull[1:]*sfull[1:], volume_profile[1:]/phi[1:])
    plt.xlabel("normalized radius")
    plt.ylabel("dV/dPsi")
    plt.title(extension)

# iota profile
iotas=vmec_data['iotas']

if do_plot_iota_profile:
    plt.figure()
    plt.plot(sfull[1:]*sfull[1:], iotas[1:], '.-')
    plt.xlabel("normalized radius")
    plt.ylabel("iotabar = iota/2pi")