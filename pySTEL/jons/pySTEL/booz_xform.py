#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 17 21:56:33 2017

@author: jonathan
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

import time

from trigfunct  import trigfunct
from transpmn import transpmn

def import_path(p):
    # check if p is already in sys.path, append if not
    try:
        sys.path.index(p)
    except ValueError:
        sys.path.append(p)

MHDpythonPath="Z:\\03_Masterarbeit\\00_programs\\MHDpython\\"
import_path(MHDpythonPath+"pyVMEC")

from read_vmec import read_vmec

verbose=False

# input parameters
mboz=35
nboz=15

extension="w7x_ref_84"
datadir="Z:\\03_Masterarbeit\\01_analysis\\24_magnetic_coordinates\\"

# construct vmec filename
vmec_filename="wout."+extension+".nc"

# read vmec data
vmec_data=read_vmec(datadir+vmec_filename)

mpol=vmec_data['mpol']
ntor=vmec_data['ntor']
ns=vmec_data['ns']
nfp=np.int(vmec_data['nfp'])
lasym=vmec_data['iasym']
mnmax=vmec_data['mnmax']
xm=vmec_data['xm']
xn=vmec_data['xn']

mpol_nyq = np.int(np.max(vmec_data['xm_nyq']))
ntor_nyq = np.int(np.max(np.abs(vmec_data['xn_nyq']))/nfp)

xm_nyq=vmec_data['xm_nyq']
xn_nyq=vmec_data['xn_nyq']
mnmax_nyq=vmec_data['mnmax_nyq']

iotas=vmec_data['iotas']

rmnc=vmec_data['rmnc']
rmns=vmec_data['rmns']
zmnc=vmec_data['zmnc']
zmns=vmec_data['zmns']
lmnc=vmec_data['lmnc']
lmns=vmec_data['lmns']
bmnc=vmec_data['bmnc']
bmns=vmec_data['bmns']

bsubumnc=vmec_data['bsubumnc']
bsubumns=vmec_data['bsubumns']
bsubvmnc=vmec_data['bsubvmnc']
bsubvmns=vmec_data['bsubvmns']

# comparison output
from netCDF4 import Dataset
filename="Z:\\03_Masterarbeit\\01_analysis\\24_magnetic_coordinates\\boozmn_w7x_ref_84.nc"
rootgrp = Dataset(filename, 'r')
booz_data=dict()

bmnc_b_import=rootgrp['/bmnc_b']
rmnc_b_import=rootgrp['/rmnc_b']
zmns_b_import=rootgrp['/zmns_b']
pmns_b_import=rootgrp['/pmns_b']
gmnc_b_import=rootgrp['/gmn_b']

#%%
# all radial surfaces
jlist=range(ns)

#     COMPUTE ACTUAL NO. THETA, PHI POINTS FOR INTEGRATIONS
#     NEEDED FOR DYNAMIC MEMORY ALLOCATION

#       js         radial point where Boozer coords. are needed
#       ns           number of vmec radial grid points
#       nu_boz       number of theta points in integration
#       nv_boz       number of zeta points in integration
#       mpol         number of theta harmonics from vmec for r,z,l
#       ntor         number of zeta harmonics from vmec (no. zeta modes = 2*ntor+1) for r,z,l
#       mpol_nyq     number of theta harmonics from vmec for bsubumn, bsubvmn
#       ntor_nyq     number of zeta harmonics from vmec (no. zeta modes = 2*ntor+1) for bsubumn, bsubvmn
#       mboz         number of boozer theta harmonics
#       nboz         number of boozer zeta harmonics
# ensure a minimul number of Boozer modes
mboz = np.max([6*mpol  , 2, mboz])
nboz = np.max([2*ntor-1, 0, nboz])
print("used # of modes for Boozer spectrum: m_boz = %d , nboz = %d"%(mboz, nboz))

nu_boz   = np.int(2*(2*mboz+1))                     #CHANGED THIS FROM 2*(3*mboz+1)
nv_boz   = np.int(2*(2*nboz+1))                     #CHANGED THIS FROM 2*(2*nboz+1)
if nboz==0: nv_boz=1
#      nu_boz   = nu_boz + MOD(nu_boz,2)          !nu_boz, nv_boz MUST be even
#      nv_boz   = nv_boz + MOD(nv_boz,2)
#nunv = nu_boz*nv_boz  # => will be overwritten later on by nu2_b*nv_boz...
mnboz = nboz+1 + (mboz-1)*(1+2*nboz)
nu2_b = np.int(nu_boz/2+1)                        #pi



# setup toridal and poloidal mode numbers
xnb=np.zeros([mnboz, 1])
xmb=np.zeros([mnboz, 1])
mnboz0 = 0
n2 = nboz
for m in range(mboz):
    n1 = -nboz
    if m==0: n1=0
    for n in range(n1, n2+1):

        if mnboz0 >= mnboz:
            print('ERROR: mnboz exceeds limit in booz_xform')

        # toroidal mode numbers
        xnb[mnboz0] = n*nfp

        # poloidal mode numbers
        xmb[mnboz0] = m

        # index counter
        mnboz0 = mnboz0 + 1
if mnboz0 != mnboz: mnboz=mnboz0


#   SCALE FACTOR FOR NORMALIZATION OF FOURIER TRANSFORMS
if lasym:
   fac = 2.0/(nu_boz*nv_boz)
else:
   fac = 2.0/((nu2_b-1)*nv_boz)

scl=np.ones([mnboz, 1])*fac
scl[np.where(np.any(np.round(xnb)==0) and np.any(np.round(xmb)==0))] = fac/2

ntorsum=np.zeros([2,1], dtype=int)
for i in range(mnmax):
    if np.round(xm[i])==0.0: ntorsum[0]+=1
    if np.round(xm[i])<=1.0: ntorsum[1]+=1


ohs = ns-1.0
hs  = 1.0/np.double(ohs)

sfull=np.sqrt(np.linspace(0,ns-1,ns)/(ns-1))
shalf=np.sqrt(np.linspace(0.5,ns-0.5,ns)/(ns-1))

#%%

#ONLY need top half of theta mesh for symmetric plasma
nu3_b=nu2_b
if lasym:
    nu3_b=nu_boz

nunv = nu3_b*nv_boz


#     COMPUTE POLOIDAL (thgrd) AND TOROIDAL (ztgrd) ANGLES
# Half-around in theta
dth=2.0*np.pi/(2.0*(nu3_b-1.0))
if lasym:
    # USE THIS FOR FULL 2-pi
    dth=2.0*np.pi/nu3_b


dzt = 2.0*np.pi/(nv_boz*nfp)

thgrd=np.zeros([nunv, 1])
ztgrd=np.zeros([nunv, 1])

lk = 0
for lt in range(nu3_b):
    for lz in range(nv_boz):
        thgrd[lk] = lt*dth
        ztgrd[lk] = lz*dzt
        lk=lk+1
        
cosm_b=np.zeros([nunv, mpol])
sinm_b=np.zeros([nunv, mpol])
cosn_b=np.zeros([nunv, ntor+1])
sinn_b=np.zeros([nunv, ntor+1])

cosm_nyq=np.zeros([nunv, mpol_nyq+1])
sinm_nyq=np.zeros([nunv, mpol_nyq+1])
cosn_nyq=np.zeros([nunv, ntor_nyq+1])
sinn_nyq=np.zeros([nunv, ntor_nyq+1])


trigfunct (thgrd, ztgrd, cosm_b,   sinm_b,   cosn_b,   sinn_b,   mpol-1,   ntor,     nunv, nfp)
trigfunct (thgrd, ztgrd, cosm_nyq, sinm_nyq, cosn_nyq, sinn_nyq, mpol_nyq, ntor_nyq, nunv, nfp)

#print("  0 <= mboz <=   %d     %d <= nboz <=   %d"%(mboz-1, -nboz, nboz))
#print("  nu_boz =    %d nv_boz =    %d"%(nu_boz, nv_boz))
#print("")
#print("             OUTBOARD (u=0)              JS          INBOARD (u=pi)")
#print("-----------------------------------------------------------------------------")
#print("  v     |B|vmec    |B|booz    Error             |B|vmec    |B|booz    Error")



pmns=np.zeros([ns, mnmax])
pmnc=np.zeros([ns, mnmax])
gpsi=np.zeros([ns, 1])
Ipsi=np.zeros([ns, 1])
jacfac=np.zeros([ns-1, 1])
    
all_tcos=np.zeros([nunv, mnmax])
all_tsin=np.zeros([nunv, mnmax])

all_cost=np.zeros([nunv, mnboz])
all_sint=np.zeros([nunv, mnboz])

cosmm=np.zeros([nunv, mboz+1])
sinmm=np.zeros([nunv, mboz+1])
cosnn=np.zeros([nunv, nboz+1])
sinnn=np.zeros([nunv, nboz+1])

start=time.time()
for mn in range(mnmax):
    
    m = np.int(np.round(xm[mn]))
    n = np.int(np.round(np.abs(xn[mn]/nfp)))
    sgn = np.sign(xn[mn])
   
    all_tcos[:, mn] = cosm_b[:,m]*cosn_b[:,n] + sinm_b[:,m]*sinn_b[:,n]*sgn
    all_tsin[:, mn] = sinm_b[:,m]*cosn_b[:,n] - cosm_b[:,m]*sinn_b[:,n]*sgn
end=time.time()
print("time for calculating all tcos, all tsin: " + np.str(end-start) + " s")

# formerly calculated using vcoords_rz
all_r12=np.zeros([ns-1, nunv])
all_z12=np.zeros([ns-1, nunv])
all_lam=np.zeros([ns-1, nunv])
all_lt =np.zeros([ns-1, nunv])
all_lz =np.zeros([ns-1, nunv])

# formerly calculated using vcoords_w
all_wp  =np.zeros([ns-1, nunv])
all_wt  =np.zeros([ns-1, nunv])
all_wz  =np.zeros([ns-1, nunv])
all_bmod=np.zeros([ns-1, nunv])

all_bbjac=np.zeros([ns-1, nunv])

uang=[]
vang=[]

#   COMPUTE BOOZER TRANSFORM, SURFACE BY SURFACE
for js in range(1,ns):
#for js in [1]:
    try:
        idx=jlist.index(js)
    except ValueError:
        next

    js=jlist[idx]

    print("\n\nworking on surface %d"%(js+1))


    # overall goal:
    #     COMPUTE FOURIER COEFFICIENTS (pmn) OF THE "SOURCE" CONTRIBUTIONS (right of Eq.10)
    #     OF THE BOOZER-TO-VMEC TRANSFORMATION FUNCTION P:
    #
    #     Theta-Booz = Theta-VMEC + Lambda + Iota*p  => eq. 3
    #     Zeta-Booz  = Zeta-VMEC  + p
    #


    # start by computing the RHS of eq. 10
    start=time.time()
    transpmn (pmns, bsubumnc, bsubvmnc, pmnc, bsubumns, bsubvmns,
          xm_nyq, xn_nyq, gpsi, Ipsi, mnmax_nyq, js, lasym)
    end=time.time()

    #
    #     BEGIN CALCULATION OF BOOZER QUANTITIES AT HALF-RADIAL
    #     MESH POINT js
    #     (ALL TRANSFORMED QUANTITIES MUST BE ON HALF-MESH FOR ACCURACY)
    #
    if (js+1 <= 1): print('js must be > 1!')
    
    scale_t1=np.zeros([mnmax, 1])    
    scale_t2=np.zeros([mnmax, 1])
    
    # even modes
    scale_t1[np.where(xm%2==0)[0]] = 1.0
    scale_t2[np.where(xm%2==0)[0]] = 1.0
      
    # odd modes
    if js>1:
        t1=1.0/sfull[js]
        t2=1.0/sfull[js-1]
    else:
        t1=1.0/sfull[1]
        t2=1.0
        rmnc[0, ntorsum[0,0]:ntorsum[1,0]] = 2.0*rmnc[1, ntorsum[0,0]:ntorsum[1,0]]/sfull[1] - rmnc[2, ntorsum[0,0]:ntorsum[1,0]]/sfull[2]
        zmns[0, ntorsum[0,0]:ntorsum[1,0]] = 2.0*zmns[1, ntorsum[0,0]:ntorsum[1,0]]/sfull[1] - zmns[2, ntorsum[0,0]:ntorsum[1,0]]/sfull[2]

    scale_t1[np.where(xm%2)[0]] = shalf[js-1]*t1
    scale_t2[np.where(xm%2)[0]] = shalf[js-1]*t2        
    
    r_t1=np.dot(all_tcos, np.multiply(rmnc[js  ,:], np.squeeze(scale_t1)))
    r_t2=np.dot(all_tcos, np.multiply(rmnc[js-1,:], np.squeeze(scale_t2)))
    
    z_t1=np.dot(all_tsin, np.multiply(zmns[js  ,:], np.squeeze(scale_t1)))
    z_t2=np.dot(all_tsin, np.multiply(zmns[js-1,:], np.squeeze(scale_t2)))
        
    # COMPUTE EVEN (in poloidal mode number) AND ODD COMPONENTS
    # OF R,Z, LAMDA IN REAL SPACE (VMEC) COORDINATES
        # on half mesh (?)
    all_r12[js-1,:] = np.mean(np.vstack([r_t1, r_t2]), axis=0)
    all_z12[js-1,:] = np.mean(np.vstack([z_t1, z_t2]), axis=0)
    
    all_lam[js-1,:] = np.dot(all_tsin, lmns[js,:])
    
    # derivates of lambda (chain rule applied to d/dPhi(lmns*sin(zeta*Phi) = lmns*zeta*cos(zeta*Phi) )
    all_lt[js-1,:]  = np.dot(all_tcos, np.multiply(lmns[js,:],      np.squeeze(xm)))
    all_lz[js-1,:]  = np.dot(all_tcos, np.multiply(lmns[js,:], -1.0*np.squeeze(xn)))
    
    
    # COMPUTE "SOURCE" PART OF TRANSFORMATION FUNCTION p==wp (RIGHT-SIDE OF EQ.10), ITS DERIVATIVES,
    # AND |B| ALL IN VMEC COORDINATES
    all_wp[js-1,:] = np.dot(all_tsin, pmns[js,:])
    
    # derivates of P
    all_wt[js-1,:]  = np.dot(all_tcos, np.multiply(pmns[js,:],      np.squeeze(xm)))
    all_wz[js-1,:]  = np.dot(all_tcos, np.multiply(pmns[js,:], -1.0*np.squeeze(xn)))
    
    all_bmod[js-1,:] = np.dot(all_tcos, bmnc[js,:])    
    
    
    
    #
    #     HERE, W IS THE PART OF THE TRANSFORMATION FUNCTION P
    #     WHICH DEPENDS ON THE SURFACE-AVERAGED COVARIANT B COMPONENTS [RIGHT-SIDE IN Eq(10)],
    #     COMPUTED IN transpmn AND vcoord_w.
    #
    #     THE FULL TRANSFORMATION FUNCTION IS THEN [Eq(10 in paper], with D=gpsi+iota*Ipsi:
    #
    #     P = (w - Ipsi*Lambda) / D
    #
    #     ALSO [Eq.(3) in paper]:
    #
    #     uboz = lambda + iota*P   (non-secular piece of boozer theta in VMEC coordinates)
    #
    #          = (gpsi*lambda + iota*w)/D
    #
    #     vboz = P                 (non-secular piece of boozer phi in VMEC coordinates)
    #
    #     FINALLY, XJAC IS THE JACOBIAN BETWEEN BOOZER, VMEC COORDINATES
    #
    #
    #     NOTE THAT LAMBDA = xl, d(lambda)/du = xlt, d(lambda)/dv = xlz
   
    jacfac[js-1] = gpsi[js] + iotas[js]*Ipsi[js]
    if (jacfac[js-1] == 0.0):
        print("ERROR: Boozer coordinate XFORM failed, jacfac = 0!")
    dem = 1.0/jacfac[js-1]
    gpsi1 = gpsi[js]*dem
    hiota1 = iotas[js]*dem
    Ipsi1 = Ipsi[js]*dem

    vboz = dem*all_wp[js-1,:] - Ipsi1*all_lam[js-1,:]       #TOTAL p in Eq.(10)
    uboz = all_lam[js-1,:] + iotas[js]*vboz
    
    psubv = dem*all_wz[js-1,:] - Ipsi1*all_lz[js-1,:] 
    psubu = dem*all_wt[js-1,:] - Ipsi1*all_lt[js-1,:]
    bsupv = 1 + all_lt[js-1,:] 
    bsupu = iotas[js] - all_lz[js-1,:]

    # Eq. (12) 
    xjac = bsupv*(1+psubv) + bsupu*psubu

#    dem = np.min(xjac)
#    dem2= np.max(xjac)
#     SAL 07/06/16
#      IF (dem*dem2 .le. zero)
#           PRINT *, ' Jacobian xjac changed sign in harfun in xbooz_xform'

#-----------------------------------------------
#     theta-boz = thgrd + uboz              
#     zeta-boz  = ztgrd + vboz
    uang = np.squeeze(thgrd)+uboz
    vang = np.squeeze(ztgrd)+vboz
    
    trigfunct (uang, vang, cosmm, sinmm, cosnn, sinnn, mboz, nboz, nunv, nfp)

    if not lasym:
        # ONLY INTEGRATE IN U HALF WAY AROUND (FOR LASYM=F)
        i = nv_boz*(nu2_b-1)                              #u=pi interval: i:imax
        imax = i+nv_boz
        for m in range(mboz+1):
            cosmm[:nv_boz,m] = 0.5*cosmm[:nv_boz,m]       #u=0
            cosmm[i:imax ,m] = 0.5*cosmm[i:imax ,m]       #u=pi
            sinmm[:nv_boz,m] = 0.5*sinmm[:nv_boz,m]       #should be zeroes
            sinmm[i:imax ,m] = 0.5*sinmm[i:imax ,m]       #should be zeroes

#     jacobian from VMEC to Boozer coords, with SPECIAL
#     radial variable s = (toroidal flux)/twopi (phip = 1)
#     cost = cos(mu-nv_boz);  sint = sin(mu-nv_boz)
    all_bbjac[js-1] = jacfac[js-1]/(all_bmod[js-1,:]*all_bmod[js-1,:])
    

#    start=time.time()
    for mn in range(mnboz):
        
        m = np.int(np.round(xmb[mn]))
        n = np.int(np.round(np.abs(xnb[mn]/nfp)))
        sgn = np.sign(xnb[mn])
       
        all_cost[:, mn] = (cosmm[:,m]*cosnn[:,n] + sinmm[:,m]*sinnn[:,n]*sgn)*xjac
        all_sint[:, mn] = (sinmm[:,m]*cosnn[:,n] - cosmm[:,m]*sinnn[:,n]*sgn)*xjac
#    end=time.time()
#    print("time for calculating all cost, all sint: " + np.str(end-start) + " s")

    bmncb = np.multiply(np.squeeze(scl), np.dot(all_bmod[js-1,:], all_cost))
    rmncb = np.multiply(np.squeeze(scl), np.dot(all_r12[js-1,:] , all_cost))
    zmnsb = np.multiply(np.squeeze(scl), np.dot(all_z12[js-1,:] , all_sint))
    pmnsb =-np.multiply(np.squeeze(scl), np.dot(vboz            , all_sint))
    gmncb = np.multiply(np.squeeze(scl), np.dot(all_bbjac[js-1] , all_cost))


    print ("diff in bmnc_boozer: " + np.str(np.max(np.abs(bmncb-bmnc_b_import[js-1]))))
    print ("diff in rmnc_boozer: " + np.str(np.max(np.abs(rmncb-rmnc_b_import[js-1]))))
    print ("diff in zmns_boozer: " + np.str(np.max(np.abs(zmnsb-zmns_b_import[js-1]))))
    print ("diff in pmns_boozer: " + np.str(np.max(np.abs(pmnsb-pmns_b_import[js-1]))))
    print ("diff in gmnc_boozer: " + np.str(np.max(np.abs(gmncb-gmnc_b_import[js-1]))))

    
    
    
    
    







