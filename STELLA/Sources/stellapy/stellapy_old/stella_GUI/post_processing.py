# -*- coding: utf-8 -*-

## Description: import variables from stella netcdf file

import sys
sys.path.append('/home/antonio/PROGRAMAS/stella/post_processing/')

from scipy.io import netcdf
import numpy as np
from scipy.integrate import simps
from stella_plots import plot_2d, movie_2d
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from PIL import *
import cv2
import jgraph

def read_stella_float(file_in,var):
    infile=file_in
    ncfile = netcdf.netcdf_file(infile,'r')   
    try:
        arr = np.copy(ncfile.variables[var][:])
        flag = True
    except KeyError:
        print('INFO: '+var+' not found in netcdf file')
        arr = np.arange(1,dtype=float)
        flag = False 
    return arr, flag


def timeavg(file_in,ft):
    infile=file_in
    ncfile = netcdf.netcdf_file(infile,'r')    
    # get time grid
    time = np.copy(ncfile.variables['t'][:])
    ntime = time.size
    it_min = ntime//3+1
    it_max = ntime-1
    it_interval = ntime - it_min
    time_steady = time[it_min:it_max]
    ntime_steady = time_steady.size
    time_interval = time[it_max-1]-time[it_min]
    favg = simps(ft[it_min:it_max],x=time_steady) \
        / time_interval
    return favg


def pp_symmetry(file_in):
    infile=file_in
    ncfile = netcdf.netcdf_file(infile,'r')
    
    # get zed grid
    zed = np.copy(ncfile.variables['zed'][:])
    nzed = zed.size
    iz0 = nzed//2+1

    # parallel velocity grid
    vpa, vpa_present = \
        read_stella_float(infile,'vpa')

    # |g|^2 averaged over kx, ky, and mu
    gzvs, gzvs_present \
        = read_stella_float(infile,'gzvs')

    cmap = 'YlGnBu'
    xlab = '$z$'
    ylab = '$v_{\parallel}$'
    title = 'avg $|g|^2$'

    gzvs_avg = np.arange(zed.size*vpa.size,dtype=float).reshape(vpa.size,zed.size)
    for i in range(gzvs.shape[1]):
        for j in range(gzvs.shape[2]):
            for k in range(gzvs.shape[3]):
                gzvs_avg[j,k] = timeavg(infile,gzvs[:,i,j,k])*np.exp(2.*vpa[j]**2)
        g_max = np.absolute(gzvs_avg).max()
        g_min = 0.0

        fig = plot_2d(gzvs_avg,zed,vpa,g_min,g_max,xlab,ylab,title+' (is= '+str(i+1)+')',cmap)

    return fig

def pp_gvmus_video(file_in,outname,ext):
    infile=file_in
    ncfile=netcdf.netcdf_file(infile,'r')

    # |g|^2 averaged over kx, ky, and z
    gvmus, gvmus_present \
        = read_stella_float(infile,'gvmus')

    # parallel velocity grid
    vpa, vpa_present = \
        read_stella_float(infile,'vpa')

    # mu grid
    mu, mu_present = \
        read_stella_float(infile,'mu')

    # get time grid
    time = np.copy(ncfile.variables['t'][:])
    ntime = time.size

    gmax = np.arange(ntime,dtype=float)
    gmin = np.arange(ntime,dtype=float)
    movie_file = outname+ '_gvmus'+ext

    for i in range(ntime):
        gmax[i] = np.absolute(gvmus[i,0,:,:].max())
    gmin[:] = 0.0
    xlabel = '$v_{\parallel}$'
    ylabel = '$\mu$'
    title = '$\int d^3 \mathbf{R} g^2$'

    movie_2d(gvmus[:,0,:,:],vpa,mu,gmin,gmax,ntime-1,movie_file,xlabel,ylabel,title,cmp='YlGnBu')
    #graph.show()


def pp_gzvs_video(file_in,outname,ext):
    infile=file_in
    ncfile=netcdf.netcdf_file(infile,'r')

    # |g|^2 averaged over kx, ky, and mu
    gzvs, gzvs_present \
        = read_stella_float(infile,'gzvs')

    # parallel velocity grid
    vpa, vpa_present = \
        read_stella_float(infile,'vpa')

    # get zed grid
    zed = np.copy(ncfile.variables['zed'][:])
    nzed = zed.size
    iz0 = nzed//2+1

    # get time grid
    time = np.copy(ncfile.variables['t'][:])
    ntime = time.size

    gmax = np.arange(ntime,dtype=float)
    gmin = np.arange(ntime,dtype=float)
    #movie_file = outname+ '_gvmus'+ext

    for i in range(ntime):
        gmax[i] = np.absolute(gzvs[i,0,:,:].max())
    gmin[:] = 0.0
    ylabel = '$v_{\parallel}$'
    xlabel = '$z$'
    title = '$\int d\mu \int d^2 \mathbf{R} g^2$'
    movie_file = outname+'_gzvs'+ext

    movie_2d(gzvs[:,0,:,:],zed,vpa,gmin,gmax,ntime-1,movie_file,xlabel,ylabel,title,cmp='YlGnBu')
    #graph.show()






















