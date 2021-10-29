import numpy as np
from stella_plots import movie_2d
from stella_data import phi2_vs_kxky, kx, ky, ntime

phi2max = np.arange(ntime,dtype=float)
phi2min = np.arange(ntime,dtype=float)
for i in range(ntime):
    phi2max[i] = np.absolute(phi2_vs_kxky[i,:,:].max())
phi2min[:] = 0.0
ylabel = '$k_x$'
xlabel = '$k_y$'
title = '$\left|\varphi(k_x,k_y)\right|^2$'
movie_file = 'phi2_vs_kxky.mp4'
movie_2d(phi2_vs_kxky,ky,kx,phi2min,phi2max,ntime-1,movie_file,xlabel,ylabel,title,cmp='YlGnBu')
