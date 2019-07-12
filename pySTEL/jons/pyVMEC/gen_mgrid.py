# -*- coding: utf-8 -*-

from pybaseutils.WS import equil_utils as equi
from pybaseutils.WS import magfield as mf
from w7x_currents import get_w7x_currents

import getpass
#%%

# reference shot
experiment_ID = "20160303.005"

# get local username
uid=getpass.getuser()

# VMECrest interface from Gavin's pybaseutils
vmec = equi.VMECrest()
vmec.verbose=False

# load currents that belong to given experiment ID
#currents = vmec.loadCoilCurrents(experiment_ID)
currents = get_w7x_currents(experiment_ID)

# generate MagneticConfiguration
magconf = mf.__makeConfigFromCurrents(currents, useIdealCoils=False, withEMload=True)

# get a new ID
mgrid_in_id = vmec.get_new_Mgrid_id('coils1152_'+experiment_ID)

# create Mgrid on server
mgrid_id = vmec.createMgrid(mgrid_in_id, magconf, minR=4.0, maxR=6.5, minZ=-0.8, maxZ=1.0, resR=125, resZ=90, resPhi=180)

# print created ID
print("created Mgrid ID \"" + mgrid_id + "\" for experiment " + experiment_ID + " for you.")