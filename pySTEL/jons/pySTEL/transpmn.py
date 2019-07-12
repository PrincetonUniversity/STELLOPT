# -*- coding: utf-8 -*-

import numpy as np
import sys

verbose=False

def transpmn(pmns, bsubumnc, bsubvmnc, pmnc, bsubumns, bsubvmns, xm, xn, gpsi, Ipsi, mnmax, js, lasym):

#     COMPUTE THE PART OF Pmns,c WHICH IS INDEPENDENT OF Lmns,c (the COVARIANT source terms
#     in Eq.10). THE NET P IN REAL SPACE IS GIVEN AS FOLLOWS:
#
#     P(FINAL) = {SUM(m,n)[pmns(m,n)*SIN(arg) + pmnc(m,n)*COS(arg)] - Ipsi*Lambda} / D
#
#     WHERE arg = m*theta - n*zeta, D = gpsi + iota*Ipsi

    for mn in range(mnmax):
        if   np.round(xm[mn]) != 0.0:
            pmns[js,mn] = bsubumnc[js,mn]/xm[mn]
            if lasym: pmnc[js,mn] = -bsubumns[js,mn]/xm[mn]
        elif np.round(xn[mn]) != 0.0:
            pmns[js,mn] = -bsubvmnc[js,mn]/xn[mn]
            if lasym: pmnc[js,mn] =  bsubvmns[js,mn]/xn[mn]
        else:
            # (m,n) = (0,0)            
            gpsi[js] = bsubvmnc[js,mn]
            Ipsi[js] = bsubumnc[js,mn]
            
            pmns[js,mn] = 0
            if lasym:
                pmnc[mn] = 0
            
#       DIAGNOSTIC: CHECK IF RADIAL CURRENT VANISHES
        #next #COMMENT THIS FOR DIAGNOSTIC DUMP
        if (verbose and (xm[mn]!=0.0 or xn[mn] != 0.0) and not (np.abs(Ipsi[js])+np.abs(gpsi[js])==0.0)):
           rad_cur = xn[mn]*bsubumnc[js,mn]+xm[mn]*bsubvmnc[js,mn]
           sys.stdout.write("js=%d, mn=%d "%(js, mn))
           rad_cur = rad_cur/(np.abs(Ipsi[js])+np.abs(gpsi[js]))
           if (np.abs(rad_cur) > 1.e-10):
               print("m: %g n: %g Radial current: %.20e"%(xm[mn],xn[mn],rad_cur))