# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 18:12:03 2017

@author: jons
"""

import numpy as np
from osa import Client
import matplotlib.pyplot as plt

# beta scan in config J index 1 (OP1.1)
# beta / % | w7x_ref_...
# 0.00     | 82
# 0.16     | 67
# 0.32     | 68
# 0.48     | 69
# 0.64     | 83
# 0.80     | 84

# URL for magnetic axis:
# http://svvmec1.ipp-hgw.mpg.de:8080/vmecrest/v1/geiger/w7x/1000_1000_1000_1000_+0390_+0390/05/0721/magneticaxis.json?phi=36

vmec = Client('http://esb.ipp-hgw.mpg.de:8280/services/vmec_v6?wsdl')

phi=36.0
def get_magaxis_radius(vmec_id):
    axis = vmec.service.getMagneticAxis(vmec_id, np.pi*phi/180.0)
    return axis.x1[0]

def get_iota_profile(vmec_id):
    iota_profile = vmec.service.getIotaProfile(vmec_id)
    return iota_profile

def get_kinetic_energy(vmec_id):
    Wkin = vmec.service.getKineticEnergy(vmec_id)
    return Wkin

    
beta=[0.00, 0.16, 0.32, 0.48, 0.64, 0.80]
magax_r = np.zeros(np.shape(beta))
Wkin = np.zeros(np.shape(beta))
iota_profiles=np.zeros([len(beta), 99])
i=0
for vmec_id in ['w7x_ref_82', 'w7x_ref_67', 'w7x_ref_68', 'w7x_ref_69', 'w7x_ref_83', 'w7x_ref_84']:
    magax_r[i] = get_magaxis_radius(vmec_id)
    Wkin[i] = get_kinetic_energy(vmec_id)
    iota_profiles[i,:]=get_iota_profile(vmec_id)
    i=i+1

#%% common plot title
plt_title="W7-X Reference equilibria conf. J index 1"

#%% Shafranov shift
plt.figure()
plt.plot(beta, magax_r, 'o-')
plt.xlabel("<beta> / %")
plt.ylabel("radial position of magnetic axis / m")
plt.title(plt_title)
plt.tight_layout()

#%% kinetic energy
plt.figure()
plt.plot(beta, Wkin/1e3, 'o-')
plt.xlabel("<beta> / %")
plt.ylabel("kinetic energy / kJ")
plt.title(plt_title)
plt.tight_layout()

#%% iota profile shift
plt.figure()
for i in np.arange(len(iota_profiles)):
    plt.plot(np.arange(0.0, 1.0, 1.0/98.0), iota_profiles[i,:], '-')
plt.xlabel("normalized radius")
plt.ylabel("iota")
plt.legend(['0.00 %', '0.16 %', '0.32 %', '0.48 %', '0.60 %', '0.80 %'], loc='lower right')
plt.title(plt_title)
plt.tight_layout()




