import numpy as np
from stella_plots import movie_2d
from stella_data import gvmus, vpa, mu, ntime

gmax = np.arange(ntime,dtype=float)
gmin = np.arange(ntime,dtype=float)
for i in range(ntime):
    gmax[i] = np.absolute(gvmus[i,0,:,:].max())
gmin[:] = 0.0
xlabel = '$v_{\parallel}$'
ylabel = '$\mu$'
title = '$\int d^3 \mathbf{R} g^2$'
movie_file = 'gvmus.mp4'
movie_2d(gvmus[:,0,:,:],vpa,mu,gmin,gmax,ntime-1,movie_file,xlabel,ylabel,title,cmp='YlGnBu')
