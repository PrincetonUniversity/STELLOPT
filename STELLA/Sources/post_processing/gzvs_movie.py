import numpy as np
from stella_plots import movie_2d
from stella_data import gzvs, vpa, zed, ntime

gmax = np.arange(ntime,dtype=float)
gmin = np.arange(ntime,dtype=float)
for i in range(ntime):
    gmax[i] = np.absolute(gzvs[i,0,:,:].max())
gmin[:] = 0.0
ylabel = '$v_{\parallel}$'
xlabel = '$z$'
title = '$\int d\mu \int d^2 \mathbf{R} g^2$'
movie_file = 'gzvs.mp4'
movie_2d(gzvs[:,0,:,:],zed,vpa,gmin,gmax,ntime-1,movie_file,xlabel,ylabel,title,cmp='YlGnBu')
