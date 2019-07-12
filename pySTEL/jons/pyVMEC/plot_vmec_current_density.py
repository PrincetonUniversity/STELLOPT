# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from read_vmec import read_vmec, finish_import, cfunct, sfunct

#woutfile = "Z:/03_Masterarbeit/data/CHECK/wout_vmec_check.nc"
woutfile = "Z:/03_Masterarbeit/data/w7x_ref_167/wout.w7x_ref_167.nc"
fingerprint_woutfile = "Z:/03_Masterarbeit/data/w7x_ref_163/wout.w7x_ref_163.nc"

vmec_data=read_vmec(woutfile)
finish_import(vmec_data)

fingerprint_vmec_data=read_vmec(fingerprint_woutfile)
finish_import(fingerprint_vmec_data)


# evaluation points
#nu = vmec_data['mpol']*4
nu=100
#nv = vmec_data['ntor']*4*vmec_data['nfp']
nv=180

theta = np.ndarray((nu,1))
zeta  = np.ndarray((nv,1))

for j in range(nu):
    theta[j]=2*np.pi*j/(nu-1)
for j in range(nv):
    zeta[j] =2*np.pi*j/(nv-1)

# flux surface geometry
r=cfunct(theta,zeta,vmec_data['rmnc'], vmec_data['xm'],vmec_data['xn'])
z=sfunct(theta,zeta,vmec_data['zmns'], vmec_data['xm'],vmec_data['xn'])

# current density on flux surfaces
sqrtg   =cfunct(theta,zeta,vmec_data['gmnc'],    vmec_data['xm_nyq'],vmec_data['xn_nyq'])
ju_sqrtg=cfunct(theta,zeta,vmec_data['currumnc'],vmec_data['xm_nyq'],vmec_data['xn_nyq'])
jv_sqrtg=cfunct(theta,zeta,vmec_data['currvmnc'],vmec_data['xm_nyq'],vmec_data['xn_nyq'])

sqrtg_fingerprint   =cfunct(theta,zeta,fingerprint_vmec_data['gmnc'],    fingerprint_vmec_data['xm_nyq'],fingerprint_vmec_data['xn_nyq'])
ju_sqrtg_fingerprint=cfunct(theta,zeta,fingerprint_vmec_data['currumnc'],fingerprint_vmec_data['xm_nyq'],fingerprint_vmec_data['xn_nyq'])
jv_sqrtg_fingerprint=cfunct(theta,zeta,fingerprint_vmec_data['currvmnc'],fingerprint_vmec_data['xm_nyq'],fingerprint_vmec_data['xn_nyq'])

#ju=ju_sqrtg/sqrtg-ju_sqrtg_fingerprint/sqrtg_fingerprint
#jv=jv_sqrtg/sqrtg-jv_sqrtg_fingerprint/sqrtg_fingerprint

ju=(ju_sqrtg-ju_sqrtg_fingerprint)/sqrtg
jv=(jv_sqrtg-jv_sqrtg_fingerprint)/sqrtg

#%%
toroidalIdx = 18
plt.figure()
plt.pcolormesh(r[:,:,toroidalIdx], z[:,:,toroidalIdx], ju[:,:,toroidalIdx])
plt.colorbar()
plt.axis("equal")
