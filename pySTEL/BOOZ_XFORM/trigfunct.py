# -*- coding: utf-8 -*-

import numpy as np

# fast sine/cosine calculation
# see p.ex. http://www.hugi.scene.org/online/coding/hugi%2016%20-%20cosine.htm

def trigfunct (uang, vang, cosm, sinm, cosn, sinn, mpol, ntor, nznt, nfp):

    cosm[:,0] = 1.0
    sinm[:,0] = 0.0
    cosm[:,1] = np.squeeze(np.cos(uang))
    sinm[:,1] = np.squeeze(np.sin(uang))
    
    for m in range(2, mpol+1):
        cosm[:,m] = cosm[:,m-1]*cosm[:,1] - sinm[:,m-1]*sinm[:,1]
        sinm[:,m] = sinm[:,m-1]*cosm[:,1] + cosm[:,m-1]*sinm[:,1]

    cosn[:,0] = 1.0
    sinn[:,0] = 0.0
    if ntor>=1:
        cosn[:,1] = np.squeeze(np.cos(vang*nfp))
        sinn[:,1] = np.squeeze(np.sin(vang*nfp))
    else:
        print("ntor < 1")

    for n in range(2, ntor+1):
        cosn[:,n] = cosn[:,n-1]*cosn[:,1] - sinn[:,n-1]*sinn[:,1]
        sinn[:,n] = sinn[:,n-1]*cosn[:,1] + cosn[:,n-1]*sinn[:,1]