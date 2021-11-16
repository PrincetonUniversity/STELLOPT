from stella_dirs import *
from stella_utils import *
from numpy import *
from pylab import *
from struct import *
from scipy import *
from scipy.special import expit
from matplotlib import *
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import LinearLocator, FixedLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from scipy.io import netcdf
import os.path
import h5py

plt.rcParams.update({'font.size': 28})
plt.rcParams['lines.linewidth'] = 2

import tabCompleter

from plotbox import *
from tabCompleter import *
import struct
import physcon as pc
import time
from os import listdir
from stella_read import *



def kspectra_movie(case, ):
    
    from stella_plots import movie_2d

    phi2_vs_kx_ky, k_x, k_y, n_time =\
                  phi2_vs_kxky(case), kx(case)[0], ky(case)[0], time(case)[1]
    
    phi2max = np.arange(n_time,dtype=float)
    phi2min = np.arange(n_time,dtype=float)
    
    for i in range(n_time):    
        phi2max[i] = np.absolute(phi2_vs_kx_ky[i,:,:].max())*rescale

    phi2min[:] = 0.0
    ylabel = '$k_x$'
    xlabel = '$k_y$'
    title  = '$|\\varphi(k_x, k_y)|^2$'

    movie_file = outdir(case) + '/phi2_vs_kxky.mp4'

    movie_2d(phi2_vs_kx_ky, k_y, k_x, phi2min,\
             phi2max, n_time-1, movie_file, xlabel, ylabel, title,cmp='YlGnBu')


def phi_at_box(case, t0=None, tube_idx=0, last=False, trange=None, merged=False, iz=None):
    
    zeta, nzed, nzed_mid = zed(case)

    k_x, nakx, nakx_mid = kx(case)
    k_y, naky = ky(case)
    
    if merged:
        time_trace, phi_data = merge(case, 'phi_vs_t')
    else:
        time_trace = time(case)[0]
        phi_data = phi_vs_t(case)
   
    n_time = size(time_trace)

    if iz == None:
        iz = nzed_mid

    phi_data = np.concatenate((phi_data[:,:,:,nakx_mid:,:,:],\
                               phi_data[:,:,:,:nakx_mid,:,:]),axis=3)
    
    if trange != None:
        phi_vs_z_time_real_to_avrg=phi_vs_z_time_real[(time_trace > trange[0]) & (time_trace < trange[1]),:]
        phi_vs_z_time_imag_to_avrg=phi_vs_z_time_imag[(time_trace > trange[0]) & (time_trace < trange[1]),:]
        phi_vs_z_real_to_plot=mean(phi_vs_z_time_real, axis=0)
        phi_vs_z_imag_to_plot=mean(phi_vs_z_time_imag, axis=0)
        
    if last == True:        
        phi_vs_z_real_to_plot = phi_vs_z_time_real[shape(phi_vs_z_time_real)[0]-1, :]
        phi_vs_z_imag_to_plot = phi_vs_z_time_imag[shape(phi_vs_z_time_imag)[0]-1, :]
        time0                 = time_trace[shape(time_trace)[0]-1]

    if t0 != None:
        phi_vs_z_real_to_cut  = phi_vs_z_time_real[(time_trace >= t0),:]
        phi_vs_z_imag_to_cut  = phi_vs_z_time_imag[(time_trace >= t0),:]
        phi_vs_z_real_to_plot = phi_vs_z_real_to_cut[0,:]
        phi_vs_z_imag_to_plot = phi_vs_z_imag_to_cut[0,:]
        time_trace            = time_trace[(time_trace >= t0)]
        time0                 = time_trace[0]




def phi_t(case, kx_idx=None, ky_idx=None, t0=None, tube_idx=0, last=False, trange=None,
          merged=True):

    zeta, nzed, nzed_mid = zed(case)

    k_x, nakx, nakx_mid = kx(case)
    k_y, naky = ky(case)
    
    if merged:
        time_trace, phi_data = merge(case, 'phi_vs_t')
    else:
        time_trace = time(case)[0]
        phi_data = phi_vs_t(case)
   
    n_time = size(time_trace)


    phi_data = np.concatenate((phi_data[:,:,:,nakx_mid:,:,:],\
                               phi_data[:,:,:,:nakx_mid,:,:]),axis=3)

    if kx_idx == None or ky_idx == None:
        print('kx or ky unspecified ==> sum over kx and ky carried out.')
        title='$\sum_{k_x,k_y}$'
        phi_vs_z_time_real = sum(sum(phi_data[:,tube_idx, :, :, :, 0],axis=2),axis=2)
        phi_vs_z_time_imag = sum(sum(phi_data[:,tube_idx, :, :, :, 1],axis=2),axis=2)

        
    elif ky_idx != None and ky_idx != None:
        print("Getting phi(z,t) for (kx, ky) = (", k_x[kx_idx], ',', k_y[ky_idx], ')')
        title = None#'t = '+str(time0) + ', (kx, ky) = (' +\
        #                str(format3(k_x[kx_idx])) + ', ' + str(format3(k_y[ky_idx])) + ')'
        # shape(phi_vs_t(case))=(N_time, N_tubes, N_zed, N_kx, N_ky, 2)
        phi_vs_z_time_real = phi_data[:,tube_idx, :, kx_idx, ky_idx, 0]
        phi_vs_z_time_imag = phi_data[:,tube_idx, :, kx_idx, ky_idx, 1]

    
    if trange != None:
        phi_vs_z_time_real_to_avrg=phi_vs_z_time_real[(time_trace > trange[0]) & (time_trace < trange[1]),:]
        phi_vs_z_time_imag_to_avrg=phi_vs_z_time_imag[(time_trace > trange[0]) & (time_trace < trange[1]),:]
        phi_vs_z_real_to_plot=mean(phi_vs_z_time_real, axis=0)
        phi_vs_z_imag_to_plot=mean(phi_vs_z_time_imag, axis=0)
        
    if last == True:        
        phi_vs_z_real_to_plot = phi_vs_z_time_real[shape(phi_vs_z_time_real)[0]-1, :]
        phi_vs_z_imag_to_plot = phi_vs_z_time_imag[shape(phi_vs_z_time_imag)[0]-1, :]
        time0                 = time_trace[shape(time_trace)[0]-1]

    if t0 != None:
        phi_vs_z_real_to_cut  = phi_vs_z_time_real[(time_trace >= t0),:]
        phi_vs_z_imag_to_cut  = phi_vs_z_time_imag[(time_trace >= t0),:]
        phi_vs_z_real_to_plot = phi_vs_z_real_to_cut[0,:]
        phi_vs_z_imag_to_plot = phi_vs_z_imag_to_cut[0,:]
        time_trace            = time_trace[(time_trace >= t0)]
        time0                 = time_trace[0]


    maxre  = phi_vs_z_real_to_plot.max()
    maxim  = phi_vs_z_imag_to_plot.max()
    maxphi = sqrt(phi_vs_z_real_to_plot**2.0+phi_vs_z_imag_to_plot**2.0).max()
        

    # Plotting |ph|
    
    ax = pl2d(xrange=[-pi,pi], yrange=[0,(maxphi**2.0)*2.0], xlabel='$\\zeta$', ylabel='$\\varphi^2$',\
         fig_size=(8.5, 7.5), title=title, ax=None, log=False)
    ax.plot(zeta, phi_vs_z_real_to_plot**2+phi_vs_z_imag_to_plot**2, '-',\
            color='black', linewidth=2, label='$\\varphi^2$')

    # Plotting Re(phi), Im(phi)
    ax = pl2d(xrange=[-pi,pi], yrange=[-maxre,maxre], xlabel='$\\zeta$', ylabel='$\\Re(\\varphi), \\Im({\\varphi})$',\
              fig_size=(8.5, 7.5), title=title,ax=None, log=False)
    ax.plot(zeta, phi_vs_z_real_to_plot, '-',\
            color='black', linewidth=2, label='$\\Re({\\varphi})$')
    ax.legend(loc=1,labelspacing=0.0, prop={'size':24})
    
    ax2 = ax.twinx()
    ax2.plot(zeta, phi_vs_z_imag_to_plot, '-', color='blue', linewidth=2, label='$\\Im(\\varphi)$')
    ax2.set_ylim([-maxim, maxim])
    ax2.legend(loc=2,labelspacing=0.0, prop={'size':24})
#    ax.autoscale()
#    ax2.autoscale()    
    show()
   

        

def phi2_kx_ky(case, trange=None, crange=[1E-3,10], log=True,\
               delta=1E-12, t0=None, last=False, movie=False, merged=False,
               cascade=True, save=False):

    from stella_plots import movie_2d

    k_x, nakx, nakx_mid = kx(case)

    if merged:
        time_trace, phi2_vs_kx_ky = merge(case, 'phi2_vs_kxky')
    else:
        time_trace = time(case)[0]
        phi2_vs_kx_ky = phi2_vs_kxky(case)
    k_y = ky(case)[0]
    
    n_time = size(time_trace)

    phi2_vs_kx_ky = np.concatenate((phi2_vs_kx_ky[:, nakx_mid:,:],\
                                    phi2_vs_kx_ky[:,:nakx_mid ,:]),axis=1)

    if trange != None:
        # phi2_vs_kx_ky has dimensions (ntime, nakx, naky)
        # Here, we remove all values of phi2_vs_kx_ky of times out of the selected
        # interval.
        phi2_vs_kx_ky_to_avrg = phi2_vs_kx_ky[(time_trace > trange[0]) & (time_trace < trange[1]),:,:]
    else:
        phi2_vs_kx_ky_to_avrg = phi2_vs_kx_ky

    if last == True:
        phi2_vs_kx_ky_avrg = phi2_vs_kx_ky_to_avrg[shape(phi2_vs_kx_ky_to_avrg)[0]-1,:,:]
    else:
        phi2_vs_kx_ky_avrg = mean(phi2_vs_kx_ky_to_avrg, axis=0)

    if t0 != None:
        phi2_vs_kx_ky_to_avrg = phi2_vs_kx_ky[(time_trace >= t0),:,:]
        time_trace            = time_trace[(time_trace >= t0)]
        phi2_vs_kx_ky_avrg    = phi2_vs_kx_ky_to_avrg[0,:,:]


    phi2_vs_ky = sum(phi2_vs_kx_ky_avrg[:,:],axis=0)
        


    ylabel = '$k_x\\rho_r$'
    xlabel = '$k_y\\rho_r$'
    title  = '$\\left<\\varphi^{2}\\right>(k_x, k_y)$ ' + 't = ' + str(time_trace[0])

    if crange != None:
        zmin = crange[0]
        zmax = crange[1]
    else:
        zmin = max(phi2_vs_kx_ky.min(), delta)
        zmax = phi2_vs_kx_ky.max()


    if save == True:
        fn = outdir(case) + '/phi2_vs_kx_ky.h5'
        hf = h5py.File(fn, 'w')
        if trange != None : hf.create_dataset('trange', data=trange)
        if last != None : hf.create_dataset('last', data=last)
        if t0 != None : hf.create_dataset('t0', data=t0)
        hf.create_dataset('kx', data=k_x)
        hf.create_dataset('ky', data=k_y)
        hf.create_dataset('phi2_vs_kx_ky', data=phi2_vs_kx_ky_avrg)
        hf.create_dataset('phi2_vs_ky', data=phi2_vs_ky)
        hf.close()
        print("File saved: ", fn)
        

    if movie:
        outfile = outdir(case) + '/phi2_vs_kxky.mp4'
        cmp='jet'
        step = 1
        fig = plt.figure(figsize=(12,8))
        x,y = np.meshgrid(k_y,k_x)

        
        im = plt.imshow(phi2_vs_kx_ky[0,:,:], cmap=cmp, vmin=zmin, vmax=zmax,
                        extent=[x.min(),x.max(),y.min(),y.max()],
                        interpolation='nearest', origin='lower', aspect='auto',\
                        norm=LogNorm())
        plt.colorbar()
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
#        plt.title(title)

        ims = []
        ims.append([im])

        for i in range(1,n_time-1,step):
            im = plt.imshow(phi2_vs_kx_ky[i,:,:], cmap=cmp, vmin=zmin, vmax=zmax,
                            extent=[x.min(),x.max(),y.min(),y.max()],
                            interpolation='nearest', origin='lower', aspect='auto',\
                            norm=LogNorm())
            #title  = '$\\left<\\varphi^{2}\\right>(k_x, k_y)$'+ 't = ' + str(time_trace[i])
#            im.title(title)
            ims.append([im])
            
        ani = animation.ArtistAnimation(fig,ims,interval=50,blit=True)
        ani.save(outfile)
        
    else:
        yran = [phi2_vs_kx_ky_avrg.min(), phi2_vs_kx_ky_avrg.max()]
        ax = cmap(xdata=k_y, ydata=k_x, zdata=phi2_vs_kx_ky_avrg,\
                  xlabel=xlabel, ylabel=ylabel, zlabel=title,\
                  ctics=False, contour=1, color=True, cmap='jet',\
                  epsname=None, crange=crange, fig_size=(8.5, 7.5), log=log)


        ax = pl2d(xrange=[k_y.min(), k_y.max()], yrange=None,\
                   xlabel='$k_y\\rho_r$', ylabel='$\\sum_{k_{x}}\left<\\varphi^{2}\\right>(k_x,k_y)$',\
                   fig_size=(8.5, 7.5), title=None, ax=None, log=False)

        xran = [k_y[1]*0.5, k_y.max()*1.5]
        yran = [0.001,20.0]
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(xran)
        ax.set_ylim(yran)

        xvec = linspace(xran[0], xran[1],50)
        yvec = xvec**(-7./3.)

        ax.plot(k_y[1:], phi2_vs_ky[1:], 's', linewidth=8, color='black', mfc='white', markersize=8)
        ax.plot(xvec, yvec, '--', color='black')
        
        show()


def get_fluxes_ky(case, trange=None, t0=None, last=False,\
                  crange=None, yrange=None, spec=0, movie=False,\
                  plot=True, save=True, verbose=True, merged=False):
    
    from stella_plots import movie_2d

    pflux_label='$\\Gamma_s/\\Gamma_{r}$' ; vflux_label='$U_{\\|,s}/U_{\\|, r}$' ; qflux_label='$Q_{s}/Q_{r}$'

    k_x, nakx, nakx_mid = kx(case)
    k_y = ky(case)[0]
    
    ns            = species(case)[1]

    if merged:
        time_trace, temp = merge(case, 'pflx_kxky')
        pflux_vs_t_kx_ky = temp[:,spec,:,:]
        dump, temp       = merge(case, 'qflx_kxky')
        qflux_vs_t_kx_ky = temp[:,spec,:,:]
        dump, temp       = merge(case, 'vflx_kxky')
        vflux_vs_t_kx_ky = temp[:,spec,:,:]
    else:
        time_trace    = time(case)[0]
        pflux_vs_t_kx_ky = pflux_vs_kxky(case)[:,spec,:,:]
        qflux_vs_t_kx_ky = qflux_vs_kxky(case)[:,spec,:,:]
        vflux_vs_t_kx_ky = vflux_vs_kxky(case)[:,spec,:,:]

    if verbose: print('Time interval = ', time_trace[0], time_trace[size(time_trace)-1
                                                                    ])
    n_time = size(time_trace)
    
    pflux_vs_t_kx_ky = np.concatenate((pflux_vs_t_kx_ky[:, nakx_mid:,:],\
                                       pflux_vs_t_kx_ky[:,:nakx_mid ,:]),axis=1)
    qflux_vs_t_kx_ky = np.concatenate((qflux_vs_t_kx_ky[:, nakx_mid:,:],\
                                       qflux_vs_t_kx_ky[:,:nakx_mid ,:]),axis=1)
    vflux_vs_t_kx_ky = np.concatenate((vflux_vs_t_kx_ky[:, nakx_mid:,:],\
                                       vflux_vs_t_kx_ky[:,:nakx_mid ,:]),axis=1)
    
    pflux_vs_t_ky    = sum(pflux_vs_t_kx_ky[:,:,:], axis=1)
    qflux_vs_t_ky    = sum(qflux_vs_t_kx_ky[:,:,:], axis=1)
    vflux_vs_t_ky    = sum(vflux_vs_t_kx_ky[:,:,:], axis=1)
    
    pflux_vs_t       = sum(pflux_vs_t_ky[:,:], axis=1)
    qflux_vs_t       = sum(qflux_vs_t_ky[:,:], axis=1)
    vflux_vs_t       = sum(vflux_vs_t_ky[:,:], axis=1)
        
    if trange != None:
        pflux_vs_t_kx_ky_to_avrg = pflux_vs_t_kx_ky[(time_trace > trange[0]) & (time_trace < trange[1]),:,:]
        pflux_vs_kx_ky_avrg      = mean(pflux_vs_t_kx_ky_to_avrg, axis=0)
        pflux_vs_kx_ky_to_plot   = pflux_vs_kx_ky_avrg

        qflux_vs_t_kx_ky_to_avrg = qflux_vs_t_kx_ky[(time_trace > trange[0]) & (time_trace < trange[1]),:,:]
        qflux_vs_kx_ky_avrg      = mean(qflux_vs_t_kx_ky_to_avrg, axis=0)
        qflux_vs_kx_ky_to_plot   = qflux_vs_kx_ky_avrg
        
        vflux_vs_t_kx_ky_to_avrg = vflux_vs_t_kx_ky[(time_trace > trange[0]) & (time_trace < trange[1]),:,:]
        vflux_vs_kx_ky_avrg      = mean(vflux_vs_t_kx_ky_to_avrg, axis=0)
        vflux_vs_kx_ky_to_plot   = vflux_vs_kx_ky_avrg
     
    elif trange == None and last == True:
        pflux_vs_kx_ky_to_plot   = pflux_vs_t_kx_ky[shape(time_trace)[0]-1,:,:]
        qflux_vs_kx_ky_to_plot   = qflux_vs_t_kx_ky[shape(time_trace)[0]-1,:,:]
        vflux_vs_kx_ky_to_plot   = vflux_vs_t_kx_ky[shape(time_trace)[0]-1,:,:]
                
    elif trange == None and t0!= None:
        pflux_vs_kx_ky_after_t0 = pflux_vs_t_kx_ky[(time_trace >= t0),:,:]
        pflux_vs_kx_ky_to_plot  = pflux_vs_kx_ky_after_t0[0,:,:]

        qflux_vs_kx_ky_after_t0 = qflux_vs_t_kx_ky[(time_trace >= t0),:,:]
        qflux_vs_kx_ky_to_plot  = qflux_vs_kx_ky_after_t0[0,:,:]

        vflux_vs_kx_ky_after_t0 = vflux_vs_t_kx_ky[(time_trace >= t0),:,:]
        vflux_vs_kx_ky_to_plot  = vflux_vs_kx_ky_after_t0[0,:,:]

    pflux_vs_ky_to_plot=sum(pflux_vs_kx_ky_to_plot, axis=0)
    qflux_vs_ky_to_plot=sum(qflux_vs_kx_ky_to_plot, axis=0)
    vflux_vs_ky_to_plot=sum(vflux_vs_kx_ky_to_plot, axis=0)
        
    if plot :
        # Plot with flux vs t
        xran = [time_trace.min(),time_trace.max()]

        yran = [pflux_vs_t.min(),pflux_vs_t.max()]       
        ax11 = pl2d(xrange=xran, yrange=yran,\
                   xlabel='$t\ (a/v_{th,r})$', ylabel=pflux_label,\
                   fig_size=(8.5, 7.5), title=None, ax=None, log=False)
        ax11.plot(time_trace, pflux_vs_t, color='black')
        
        yran = [qflux_vs_t.min(),qflux_vs_t.max()]
        ax12 = pl2d(xrange=xran, yrange=yran,\
                   xlabel='$t\ (a/v_{th,r})$', ylabel=qflux_label,\
                   fig_size=(8.5, 7.5), title=None, ax=None, log=False)
        ax12.plot(time_trace, qflux_vs_t, color='black')

        yran = [vflux_vs_t.min(),vflux_vs_t.max()]
        ax13 = pl2d(xrange=xran, yrange=yran,\
                   xlabel='$t\ (a/v_{th,r})$', ylabel=vflux_label,\
                   fig_size=(8.5, 7.5), title=None, ax=None, log=False)
        ax13.plot(time_trace, vflux_vs_t, color='black')


        if trange!= None:
            ax11.barh(y=yran[0], width=(trange[1]-trange[0]),\
                     height=yran[1]-yran[0], \
                     left=trange[0], align='edge', facecolor='yellow', alpha=0.3)
            ax12.barh(y=yran[0], width=(trange[1]-trange[0]),\
                     height=yran[1]-yran[0], \
                     left=trange[0], align='edge', facecolor='yellow', alpha=0.3)
            ax13.barh(y=yran[0], width=(trange[1]-trange[0]),\
                      height=yran[1]-yran[0], \
                      left=trange[0], align='edge', facecolor='yellow', alpha=0.3)
                        
        # Plot with flux vs kx ky
        if crange == None:
            c1range = [pflux_vs_kx_ky_to_plot.min(),pflux_vs_kx_ky_to_plot.max()]
            c2range = [qflux_vs_kx_ky_to_plot.min(),qflux_vs_kx_ky_to_plot.max()]
            c3range = [vflux_vs_kx_ky_to_plot.min(),vflux_vs_kx_ky_to_plot.max()]
        else:
            c1range = crange
            c2range = crange
            c3range = crange
        
        ax21 = cmap(k_y, k_x, pflux_vs_kx_ky_to_plot[:,:],\
                    xlabel='$k_y\\rho_{r}$', ylabel='$k_x\\rho_{r}$',zlabel=pflux_label,\
                    crange=c1range, log=False)
        ax22 = cmap(k_y, k_x, qflux_vs_kx_ky_to_plot[:,:],\
                    xlabel='$k_y\\rho_{r}$', ylabel='$k_x\\rho_{r}$',zlabel=qflux_label,\
                    crange=c2range, log=False)
        ax23 = cmap(k_y, k_x, vflux_vs_kx_ky_to_plot[:,:],\
                    xlabel='$k_y\\rho_{r}$', ylabel='$k_x\\rho_{r}$',zlabel=vflux_label,\
                    crange=c3range, log=False)
                
        # Plot with flux vs ky
        ax31 = pl2d(xrange=[k_y.min(),k_y.max()],\
                    yrange=[pflux_vs_ky_to_plot.min(), pflux_vs_ky_to_plot.max()],\
                    xlabel='$k_{y}\\rho_{r}$', ylabel=pflux_label,\
                    fig_size=(8.5, 7.5), title=None, ax=None, log=False)
        
        ax32 = pl2d(xrange=[k_y.min(),k_y.max()],\
                    yrange=[qflux_vs_ky_to_plot.min(), qflux_vs_ky_to_plot.max()],\
                    xlabel='$k_{y}\\rho_{r}$', ylabel=qflux_label,\
                    fig_size=(8.5, 7.5), title=None, ax=None, log=False)
        
        ax33 = pl2d(xrange=[k_y.min(),k_y.max()],\
                    yrange=[vflux_vs_ky_to_plot.min(), vflux_vs_ky_to_plot.max()],\
                    xlabel='$k_{y}\\rho_{r}$', ylabel=vflux_label,\
                    fig_size=(8.5, 7.5), title=None, ax=None, log=False)
        
        ax31.plot(k_y, pflux_vs_ky_to_plot, color='black')
        ax32.plot(k_y, qflux_vs_ky_to_plot, color='red')
        ax33.plot(k_y, vflux_vs_ky_to_plot, color='purple')
                
        show()

    if save == True:
        fn = outdir(case) + '/fluxes_spec_'+str(spec)+'_kx_ky.h5'
        hf = h5py.File(fn, 'w')
        if trange != None : hf.create_dataset('trange', data=trange)
        if last != None : hf.create_dataset('last', data=last)
        hf.create_dataset('spec', data=spec)
        if t0 != None : hf.create_dataset('t0', data=t0)
        hf.create_dataset('kx', data=k_x)
        hf.create_dataset('ky', data=k_y)
        hf.create_dataset('pflux_vs_kx_ky', data=pflux_vs_kx_ky_to_plot)
        hf.create_dataset('pflux_vs_ky',    data=pflux_vs_ky_to_plot)
        hf.create_dataset('qflux_vs_kx_ky', data=qflux_vs_kx_ky_to_plot)
        hf.create_dataset('qflux_vs_ky',    data=qflux_vs_ky_to_plot)
        hf.create_dataset('vflux_vs_kx_ky', data=vflux_vs_kx_ky_to_plot)
        hf.create_dataset('vflux_vs_ky',    data=vflux_vs_ky_to_plot)
        hf.close()
        print("File saved: ", fn)

        

#def merge_fluxes(case):
    # This function localizes recursively all the *fluxes files within
    # a case directory.

def density_kx_ky(case, trange=None, crange=[1E-3,10], log=True,\
               delta=1E-12, t0=None, last=False, movie=False, merged=False,
               cascade=True, save=False, spec=0, tube=0, zpos=0):

    if merged:
        time_trace, density = merge(case, 'density')
    else:
        time_trace = time(case)[0]
        density = density_vs_kxky(case)
    
    dn_kx_ky = density[:,spec,tube,zpos,:,:,0]+\
               density[:,spec,tube,zpos,:,:,1]*1j
    dn2_kx_ky = abs(dn_kx_ky)**2.0
    
    k_x, nakx, nakx_mid = kx(case)
    k_y = ky(case)[0]
    n_time = size(time_trace)

    dn2_vs_kx_ky = np.concatenate((dn2_kx_ky[:, nakx_mid:,:],\
                                   dn2_kx_ky[:,:nakx_mid ,:]),axis=1)

    if trange != None:
        # dn2_vs_kx_ky has dimensions (ntime, nakx, naky)
        # Here, we remove all values of dn2_vs_kx_ky of times out of the selected
        # interval.
        dn2_vs_kx_ky_to_avrg = dn2_vs_kx_ky[(time_trace > trange[0]) & (time_trace < trange[1]),:,:]
    else:
        dn2_vs_kx_ky_to_avrg = dn2_vs_kx_ky

    if last == True:
        dn2_vs_kx_ky_avrg = dn2_vs_kx_ky_to_avrg[shape(dn2_vs_kx_ky_to_avrg)[0]-1,:,:]
    else:
        dn2_vs_kx_ky_avrg = mean(dn2_vs_kx_ky_to_avrg, axis=0)

    if t0 != None:
        dn2_vs_kx_ky_to_avrg = dn2_vs_kx_ky[(time_trace >= t0),:,:]
        time_trace           = time_trace[(time_trace >= t0)]
        dn2_vs_kx_ky_avrg    = dn2_vs_kx_ky_to_avrg[0,:,:]


    dn2_vs_ky = sum(dn2_vs_kx_ky_avrg[:,:],axis=0)
        


    ylabel = '$k_x\\rho_r$'
    xlabel = '$k_y\\rho_r$'
    title  = '$\\left|delta n(k_x, k_y)\\right|^2$'

    if crange != None:
        zmin = crange[0]
        zmax = crange[1]
    else:
        zmin = max(dn2_vs_kx_ky.min(), delta)
        zmax = dn2_vs_kx_ky.max()


    if save == True:
        fn = outdir(case) + '/dn2_vs_kx_ky.h5'
        hf = h5py.File(fn, 'w')
        if trange != None : hf.create_dataset('trange', data=trange)
        if last != None : hf.create_dataset('last', data=last)
        if t0 != None : hf.create_dataset('t0', data=t0)
        hf.create_dataset('kx', data=k_x)
        hf.create_dataset('ky', data=k_y)
        hf.create_dataset('dn2_vs_kx_ky', data=dn2_vs_kx_ky_avrg)
        hf.create_dataset('dn2_vs_ky', data=dn2_vs_ky)
        hf.close()
        print("File saved: ", fn)
        

    if plot:
        yran = [dn2_vs_kx_ky_avrg.min(), dn2_vs_kx_ky_avrg.max()]
        ax = cmap(xdata=k_y, ydata=k_x, zdata=dn2_vs_kx_ky_avrg,\
                  xlabel=xlabel, ylabel=ylabel, zlabel=title,\
                  ctics=False, contour=1, color=True, cmap='jet',\
                  epsname=None, crange=crange, fig_size=(8.5, 7.5), log=log)


        ax = pl2d(xrange=[k_y.min(), k_y.max()], yrange=None,\
                   xlabel='$k_y\\rho_r$', ylabel='$\\sum_{k_{x}}\left<\\delta n^{2}\\right>(k_x,k_y)$',\
                   fig_size=(8.5, 7.5), title=None, ax=None, log=False)

        xran = [k_y[1]*0.5, k_y.max()*1.5]
        yran = [0.001,20.0]
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(xran)
        ax.set_ylim(yran)

        xvec = linspace(xran[0], xran[1],50)
        yvec = xvec**(-7./3.)

        ax.plot(k_y[1:], dn2_vs_ky[1:], 's', linewidth=8, color='black', mfc='white', markersize=8)
        ax.plot(xvec, yvec, '--', color='black')
        
        show()



    
def get_density(case, trange=None, tposition=None, zposition=0,\
                t0=None, last=False, crange=None, yrange=None, spec=0,\
                movie=True, tube=0, ikx=None, iky=None, plot=True):
    #
    # ('density', <class 'netCDF4._netCDF4.Variable'>
    #  float64 density(t, species, tube, zed, theta0, ky, ri)
    #  long_name: perturbed density vs (ky,kx,z,t)
    #  unlimited dimensions: t
    #  current shape = (ntime,nspec, ntubes, nzed, nakx, naky, 2)
    #
    #  Note that: for even nzed, the code considers the next odd integer
    #             for the dimention of zed.

    Lx = lx(case)
    Ly = ly(case)
    kx_vec=kx_stella(case)
    ky_vec=ky(case)[0]    
    kx_min, ky_min = kx_vec[1], ky_vec[1]
    
    Nkx, Nky  = size(kx_vec), size(ky_vec)
    
    x = linspace(-Lx/2, Lx/2, Nkx)
    y = linspace(-Ly/2, Ly/2, Nky)   
    X, Y = meshgrid(x, y)

    dn_x_y = empty((Nky,Nkx), dtype='cfloat')

    f, ax = plt.subplots(figsize=(8.5, 7.5))
    ax.set_ylabel('$x$ [$\\rho_\\mathrm{i}$]')
    ax.set_xlabel('$y$ [$\\rho_\\mathrm{i}$]')
    ax.axis([x[0], x[size(x)-1], y[0], y[size(y)-1]])
    
    for i in arange(0,time(case)[1]):

        dn_kx_ky = density_vs_kxky(case)[i,spec,tube,zposition,:,:,0]+\
                   density_vs_kxky(case)[i,spec,tube,zposition,:,:,1]*1j

        dn_x_y=0

        for i in arange(Nkx):
            for j in arange(Nky):
                dn_x_y = dn_x_y+dn_kx_ky[i,j]*exp(1j*kx_vec[i]*X+1j*ky_vec[j]*Y)
                
        ax.clear()
#        cbar.clear()
        ax.pcolor(x, y,real(dn_x_y), cmap='seismic', shading='auto')
#        cbar = plt.colorbar(img) 
#        cbar.set_label('$\\delta n$')
        plt.pause(0.01)
        draw()

##    #
##    # Thus, we have to transpose it as ifft2
##    # gets as input F(S)(Nky, Nkx)
##    #
##    density_x_y   = ifft2(density_kx_ky.T)
##    #
##    # print(shape(density_x_y)) ==> (Ny=Nky, Nx=Nkx), OK!
##    X, Y = meshgrid(linspace(-pi, pi), linspace(-2*pi, 2*pi))

##    fig, ax1 = plt.subplots(figsize=(12, 6))
##    pos = ax1.imshow(real(density_x_y), interpolation=None, cmap=cm.RdYlGn,
##                     origin='lower', extent=[-Lx/2, Lx/2, -Ly/2, Ly/2])
##    fig.colorbar(pos, ax=ax1, cmap='seismic')
    show()



def fft2prueba():
    X, Y = meshgrid(linspace(-pi,pi,20), linspace(-2*pi,2*pi,30))

    F = cos(2*X)+1

    imshow(imag(fft2(F)))

#    imshow(real(ifft2(fft2(F))))

    show()


def kx_min(case):
    kx_mid = kx(case)[2]
    return kx(case)[0][kx_mid]

def ky_min(case):
    return ky(case)[0][1]

def ly(case):
    return 2*pi/ky_min(case)

def lx(case):
    return 2*pi/kx_min(case)

def get_moment(case, trange=None, t0=None, last=False,\
               crange=None, yrange=None, moment='n', spec=0, movie=True,\
               tube=0, ikx=None, iky=None, plot=True):

    # Within *out.nc we have: density(t, nspec, ntubes, nztot, nakx, naky, 2==>complex)
    #                         upar(t, nspec, ntubes, nztot, nakx, naky, 2==>complex)
    #                         temperature(t, nspec, ntubes, nztot, nakx, naky, 2==>complex)
    from stella_plots import movie_2d

    if moment == 'n': moment_label='$\\delta{n}^2_s$'
    if moment == 'u': moment_label='$\\delta{U}^2_{\\|,s}$'
    if moment == 't': moment_label='$\\delta{T}^2_{s}$'

    k_x, nakx, nakx_mid = kx(case)
    k_y, n_time   = ky(case)[0], time(case)[1]
    time_trace    = time(case)[0]
    zeta          = zed(case)[0]

    # We select the species and the tube ==> (t,nztot,nakx,naky, 2)
    if moment == 'n': moment_vs_t_zeta_kx_ky = density_vs_kxky(case)[:,spec,tube, :,:,:]
    if moment == 'u': moment_vs_t_zeta_kx_ky = upar_vs_kxky(case)[:,spec,tube, :,:,:]
    if moment == 't': moment_vs_t_zeta_kx_ky = temperature_vs_kxky(case)[:,spec,tube,:,:,:]
    
    moment_vs_t_zeta_kx_ky = np.concatenate((moment_vs_t_zeta_kx_ky[:, nakx_mid:,:,:],\
                                             moment_vs_t_zeta_kx_ky[:,:nakx_mid ,:,:]),axis=1)
    if ikx == None and iky == None:
        moment2_vs_t_zeta_kx_ky = moment_vs_t_zeta_kx_ky[:,:,:,:,0]**2.0+\
                                  moment_vs_t_zeta_kx_ky[:,:,:,:,1]**2.0
        moment2_vs_t_zeta_ky    = sum(moment2_vs_t_zeta_kx_ky[:,:,:,:], axis=2)
        moment2_vs_t_zeta       = sum(moment2_vs_t_zeta_ky[:,:,:], axis=2)
        mean_moment2_vs_zeta    = mean(moment2_vs_t_zeta[:,:], axis=0)

    elif ikx != None and iky != None:
        print('Selecting the mode with kx, ky =', k_x[ikx], k_y[iky])
        moment2_vs_t_zeta  = moment_vs_t_zeta_kx_ky[:,:,ikx,iky,0]#**2.0+\
        #                             moment_vs_t_zeta_kx_ky[:,:,ikx,iky,1]**2.0
        moment2_vs_t_zeta  = moment_vs_t_zeta_kx_ky[:,:,ikx,iky,1]
        
        mean_moment2_vs_zeta    = mean(moment2_vs_t_zeta[:,:], axis=0)        

    if trange != None:
        moment2_vs_t_zeta_to_avrg = moment2_vs_t_zeta[(time_trace > trange[0]) & (time_trace < trange[1]),:]#-\
#                                    mean_moment2_vs_zeta
        moment2_vs_zeta_avrg      = mean(moment2_vs_t_zeta_to_avrg, axis=0)
        moment2_vs_zeta_to_plot   = moment2_vs_zeta_avrg
        
    elif trange == None and last == True:
        moment2_vs_zeta_to_plot = moment2_vs_t_zeta[shape(time_trace)[0]-1,:]
        
    elif trange == None and t0!= None:
        moment2_vs_zeta_after_t0 = moment2_vs_t_zeta[(time_trace >= t0),:]
        moment2_vs_zeta_to_plot  = moment2_vs_zeta_after_t0[0,:]

    if plot:
        xran = [zeta.min(), zeta.max()]
        yran = [moment2_vs_zeta_to_plot.min(), moment2_vs_zeta_to_plot.max()]

        ax1 = pl2d(xrange=xran, yrange=None, xlabel='$z$', ylabel=moment_label,\
                   fig_size=(8.5, 7.5), title=None, ax=None, log=False)
        ax1.autoscale()
        ax1.plot(zeta, moment2_vs_zeta_to_plot, '-', color='navy',\
                 linewidth=4, markerfacecolor='white')
        show()
        
    return moment_vs_t_zeta_kx_ky

def get_kperp2(case):
    # This function builds kperp2 as done internally by the code
    # See dist_fn.f90
    #
    # do iky = 1, naky
    #   if (zonal_mode(iky)) then
    #      do ikx = 1, nakx
    #         kperp2(iky,ikx,:,:) = akx(ikx)*akx(ikx)*gds22/(geo_surf%shat**2)
    #      end do
    #   else
    #      do ikx = 1, nakx
    #         kperp2(iky,ikx,:,:) = aky(iky)*aky(iky) &
    #              *(gds2 + 2.0*theta0(iky,ikx)*gds21 &
    #              + theta0(iky,ikx)*theta0(iky,ikx)*gds22)
    #      end do
    #   end if
    # end do
    kx                = kx_stella(case)
    ky                = ky_stella(case)
    z, nzed, mid_nzed = zed(case)
    
    kperp2 = empty((size))



def omega_k(case, last=False, view=True, write=False):
    # case should be a list of input files.
    for i in case:
        n_ky, n_kx = ky(i)[1], n_kx = kx(i)[1]

        if n_ky == 1 and n_kx == 1:
            omega_file = i+'.omega'
            omega_data = loadtxt(file_omega, dtype='float')
            omega_data_finite = omega_data.dropna()
#            omega_data_last   = omega_data_finite[]
            
            

def omega(case, last=False, view=True, yrange=None, xrange=None):
    # Note that:
    # time   ky   kx   Re[om]   Im[om]   Re[omavg]  Im[omavg]
    
    for i in arange(0,size(case)):
        datafile = outfile(case[i],'omega')
        n_time   = time(case[i])[1] # Including t=0
        n_ky     = ky(case[i])[1]
        n_kx     = kx(case[i])[1]

        data     = loadtxt(datafile, dtype='float').reshape(n_time-1, n_kx, n_ky, 7)

        if view:
            pl2y(ky(case[i])[0], y1data=data[n_time - 2, 0,:,3], y2data=data[n_time - 2, 0,:,4],\
                 xlabel='$k_y\\rho_{i}$', ylabel='$\\omega a/v_{\mathrm{th},\mathrm{i}}$',\
                 key1='$a\\Re(\\omega)/v_{\mathrm{th},\mathrm{i}}$',\
                 key2='$a\\Im(\\omega)/v_{\mathrm{th},\mathrm{i}}$', xrange=[0,20],\
                 yrange=yrange,fig_size=(8.5, 7.5), wp=1, ax=None, mkt="o", mkc='red',\
                 title="$k_x="+str(kx(case[i])[0][0])+"$",\
                 hline1=None, hline2=None, vshadow=None, ls1=1, ls2=2)
            show()
            
        if n_kx == 1 and n_ky > 1:
                # The evolution of omega(ky, t) can be represented
                t      = time(case[i])[0][1:]
                kyrhoi = ky(case[i])[0]
                ky_re  = data[:, 0,:,3]
                ky_im  = data[:, 0,:,4]
                [X,Y]  = meshgrid(t, kyrhoi)
                Z      = ky_im.T

                surf(X, Y, Z, xlabel='$t$', ylabel='$k_y\\rho_{i}$',\
                     zlabel='$a\\Im(\\omega)/v_{\mathrm{th},\mathrm{i}}$',\
                     title="$k_x="+str(kx(case[i])[0][0])+"$", zrange=yrange, ax=None)
                show()

        if n_kx == 1 and n_ky == 1:
                # The evolution of omega(t) can be represented
                t      = time(case[i])[0][1:]
                kyrhoi = ky(case[i])[0]
                om_re  = data[:, 0,0,3]
                notnan = logical_not(isnan(om_re))
                om_re  = om_re[notnan]
                om_im  = data[:, 0,0,4]
                om_im  = om_im[notnan]
                t      = t[notnan]
                
                xlims = [t.min(), t.max()]
                ylims = [min(om_im.min(), abs(om_re).min())*1.05, max(om_im.max(), abs(om_re).max())*1.05]

                if view: 
                    pl2y(t, abs(om_re), om_im, xlabel='$t$', ylabel='$\\omega a/v_{\mathrm{th},\mathrm{i}}$',\
                         title="$(k_x, k_y)=("+str(kx(case[i])[0][0])+", " + str(ky(case[i])[0][0])+")$",\
                         xrange=xlims, yrange=ylims,\
                         key1='$a\\Re(\\omega)/v_{\mathrm{th},\mathrm{i}}$',\
                         key2='$a\\Im(\\omega)/v_{\mathrm{th},\mathrm{i}}$',\
                         ax=None)

                

                    show()
            
    if last:
        return data[n_time-2, :, :, :]
    else:
        return data

def omega_all(case, last=False, view=True, yrange=None, xrange=None, scan='ky'):
    # This function takes a folder where various single (kx,ky)
    # runs have been performed and gather all results in one file,
    # taking the last instant with finite values of omega.
#    hidden = casestr(case).startswith('.')

    input_all=[outdir(case) + '/' + i for i in inputlist(case)]

    for i in input_all:
        # We assume that if no omega file related to a certain input is present
        # then the case is not valid.
        #print(outfile(case + '/' + i, quant='omega'))
        if not os.path.isfile(outfile(i, quant='omega')):
            print("No *omega data found for case with input " + i +'. Case ignored.') 
            input_all.remove(i)
            
    nval            = size(input_all)
    ns              = nspec([i])
    omega_all_data  = empty((nval, 7))
    phi2_all_data   = empty(nval)

    
    for i in arange(0, nval-1):
        omega_file = outfile(input_all[i], quant='omega')
        omega_all_data[i,:]=omega_last(omega_file)
        phi2_all_data[i]=phi2_vs_kxky(input_all[i])[shape(phi2_vs_kxky(input_all[i]))[0]-1,0,0]
        print(phi2_all_data[i])

    
    

def omega_last(omegafile):
    # This function returns the last row with finite values of omega
    # of a certaing omega file. The full path must be pass onto this function.
    # It assumes that only pair (kx,ky) are considered
    omegadata   = loadtxt(omegafile, dtype='float')
    notnan      = isfinite(omegadata[:,3])
    omegadata_f = omegadata[notnan,:]
    return omegadata_f[shape(omegadata_f)[0]-1,:]
        
        


def multi_lin_data(case=None, view=True, yrange=None, xrange=None, sort='ky',\
                   off=[], save=True, ql='per', last=True):
    # Note that:
    # ql == 'mik' applies formula (1) of Mikkelen et al. PoP 21 082302 (2014)
    # ql == 'per' applies formula (13)-(14) of Helander and Zocco PPCF 60 084006 (2018)
    
    # time   ky   kx   Re[om]   Im[om]   Re[omavg]  Im[omavg]

    if size(case) == 1:
        # When a directory is pass onto "case", several runs of single kx,ky
        # are assumed. The input files on that directory are read, neglecting
        # those with no *out.nc output.
        cases=[case[0] + '/' + i for i in inputlist(case[0])]
    else:
        # When a list of input files of interest is pass onto, we take
        # it as a set of cases to run the diagnostics ql_fluxes over.
        cases = case

    # We count first the number of input files with existing *out.nc
    count=0
    not_valid=[]

    print("\nChecking existence of *out.nc file corresponding inputs.")

    for i in cases:
        # Analysis is performed only if *out.nc file is present
        if os.path.isfile(outfile(i, 'out.nc')):           
            count = count + 1
            ns    = nspec(i)
        else:
            print("Case ", i, " removed from input list.")
            not_valid.append(i)
            
        for j in off:
            if j in i: not_valid.append(i)
            
    print("Number of found input files with corresponding output *.nc =", count)
            
    for j in not_valid:
        cases.remove(j)

    # Columns will contain kxrhoi, kyrhoi, re(omega) im(omega) nspec*3 fluxes
    ncol = 7 + ns*3
    lin_data = empty((size(cases), ncol), dtype='float')
    
    zeta_norm, nzed, iz0 = zed(cases[0])
    phi2data = empty((size(cases), nzed), dtype='float')


    for i in cases:
        # Analysis is performed only if *out.nc file is present
        print("ql_fluxes diagnostics running over case", i)
        lin_data[cases.index(i),:] = ql_fluxes(case=i, view=False, formula=ql)[1]
        ntime                      = shape(phi_vs_t(i))[0]
        phi2data[cases.index(i),:]  = phi_vs_t(i)[ntime-1,0,:,0,0,0]**2.0+\
                                      phi_vs_t(i)[ntime-1,0,:,0,0,1]**2.0
        phi2data[cases.index(i),:]  = phi2data[cases.index(i),:]/phi2data[cases.index(i),:].max()
         
        if lin_data[cases.index(i),6] < 100:
            print("Warning! phi2 found too low: setting omega, gamma and fluxes to NaN.")
            lin_data[cases.index(i),3]  = NaN
            lin_data[cases.index(i),2]  = NaN
            lin_data[cases.index(i),7:] = NaN
            phidata = NaN


    # We sort the data by increasing ky order
    olab   = '$\\omega a/v_{\\mathrm{th},i}$'
    glab   = '$\\gamma a/v_{\\mathrm{th},i}$'
    qlGlab = '$\\Gamma^{\\mathrm{ql}}$ [a.u.]'
    qlUlab = '$U^{\\mathrm{ql}}$ [a.u.]'
    qlQlab = '$Q^{\\mathrm{ql}}$ [a.u.]'
   
    
    if sort == 'ky':
        col = 1 ; xlab = '$k_{y}\\rho_i$'
    if sort == 'kx':
        col = 0 ; xlab = '$k_{x}\\rho_i$'

    sort     = lin_data[:,col].argsort()
    lin_data = lin_data[sort]
    phi2data = phi2data[sort]

    if save == True:
        fn = outdir(cases[0]) + '/multi_lin_data_'+ql+'.h5'
        hf = h5py.File(fn, 'w')
        hf.create_dataset('kx', data=lin_data[:,0])
        hf.create_dataset('ky', data=lin_data[:,1])
        hf.create_dataset('omega_ql_last', data=lin_data[:,2])
        hf.create_dataset('gamma_ql_last', data=lin_data[:,3])
        hf.create_dataset('omega_full_last', data=lin_data[:,4])
        hf.create_dataset('gamma_full_last', data=lin_data[:,5])        
        hf.create_dataset('phi2_ql_last', data=lin_data[:,6])
        hf.create_dataset('qlpflx', data=lin_data[:,7:7+ns])
        hf.create_dataset('qlvflx', data=lin_data[:,7+ns:7+2*ns])
        hf.create_dataset('qlqflx', data=lin_data[:,7+2*ns:7+3*ns])
        hf.create_dataset('vth_r', data=lin_data[:,7+2*ns:7+3*ns])
        hf.create_dataset('phi2_vs_ky', data=phi2data)
        hf.close()
        print("File saved: ", fn)

    if view == True:

        # Plotting
        xran = [min(lin_data[:,col]), max(lin_data[:,col])]

        # omega
        ax1 = pl2d(xrange=xran, yrange=None, xlabel=xlab, ylabel=olab,\
                   fig_size=(8.5, 7.5), title=None, ax=None, log=False)
        ax1.autoscale()
        ax1.plot(lin_data[:,col], lin_data[:,2], 'o-', color='navy',\
                 linewidth=4, markerfacecolor='white')
        if last == True:
            ax1.plot(lin_data[:,col], lin_data[:,4], 'o-', color='gray',\
                     linewidth=2, markerfacecolor='white')
        # gamma
        ax2 = pl2d(xrange=xran, yrange=None, xlabel=xlab, ylabel=glab,\
                   fig_size=(8.5, 7.5), title=None, ax=None, log=False)
        ax2.autoscale()
        ax2.plot(lin_data[:,col], lin_data[:,3], 'o-', color='crimson',\
                 linewidth=4, markerfacecolor='white')
        if last == True:
            ax2.plot(lin_data[:,col], lin_data[:,5], 'o-', color='gray',\
                     linewidth=2, markerfacecolor='white')            

        
        for k in arange(0, ns):
            fig = plt.figure(figsize=(18, 9))
            fig.subplots_adjust(left=0.1, wspace=0.35)

            title = 'species '+str(k+1)
            ax1 = fig.add_subplot(131)
            pl2d(xrange=xran, yrange=None, xlabel=xlab, ylabel=qlGlab,\
                 fig_size=(8.5, 7.5), title=title, ax=ax1)
            ax1.autoscale()
            ax1.plot(lin_data[:,col], lin_data[:,7+k], 'o-', color='darkviolet',\
                     linewidth=4, markerfacecolor='white')
            
            ax2 = fig.add_subplot(132)
            pl2d(xrange=xran, yrange=None, xlabel=xlab, ylabel=qlUlab,\
                 fig_size=(8.5, 7.5), title=title, ax=ax2)
            ax2.autoscale()
            ax2.plot(lin_data[:,col], lin_data[:,7+ns+k], 'o-', color='seagreen',\
                     linewidth=4, markerfacecolor='white')

            ax3 = fig.add_subplot(133)
            pl2d(xrange=xran, yrange=None, xlabel=xlab, ylabel=qlQlab,\
                 fig_size=(8.5, 7.5), title=title, ax=ax3)
            ax3.autoscale()
            ax3.plot(lin_data[:,col], lin_data[:,7+2*ns+k], 'o-', color='darkorange',\
                     linewidth=4, markerfacecolor='white')

        
        cmap(xdata=zeta_norm, ydata=lin_data[:,col], zdata=phi2data,\
             xlabel='$z$', ylabel=xlab, zlabel='$\\varphi^2/\\varphi^{2}_{\\mathrm{max}}$',\
             ctics=False, num=None, contour=1, color=True, cmap='jet', vfield=False, vxdata=None,\
             vydata=None, cont2=None, epsname=None, crange=None, fig_size=(8.5, 7.5),\
             log=False, ax=None)
            
        show()


#            print(type(ks), type(fluxes))
            

                                 

#            multi_lin_data  = concatenate([array(fluxes[0],fluxes[1]), fluxes[2]])
#            print(multi_lin_data)
#            omega  = 

def fluxes(case=None, plot=False, si=False, trange=None, tref=None, sign=1, nref=None,\
           merged=False, spec=0, save=True):
    i        = case
    fluxfile = outfile(i,'fluxes')
    
    if (merged):
        fluxdata = merge(case, quant='fluxes')[1]
        #        fluxfile = fluxfile + '_merged'

    else:
        fluxdata = loadtxt(fluxfile, dtype='float')
        
    nspecies = species(i)[1]
    
    if si:
        time_ref      = 1/ref_values(case, tref=tref)[5]
        fluxes_ref    = array(ref_values(case, tref=tref, nref=nref)[9:12])

        # We undo the normalization to get SI units
        fluxdata[:,0]=fluxdata[:,0]*time_ref
        for j in arange(0,3):
            #            for k in arange(0,nspecies):
            #fluxdata[:,1+j*3+k]=sign*fluxdata[:,1+j*3+k]*fluxes_ref[j]
            fluxdata[:,1+j*nspecies+spec]=sign*fluxdata[:,1+j*3+spec]*fluxes_ref[j]
            #fluxdata[:,1+spec*3+j]=sign*fluxdata[:,1+spec*3+j]*fluxes_ref[j]
            
    if trange :
        flux_data_to_avrg     = fluxdata[(fluxdata[:,0]>trange[0]) & (fluxdata[:,0]< trange[1]) ,:]
        
#        flux_data[(flux_data[:,0]>t0) & (flux_data[:,0]< tf) ,:]
        fluxes_avrg = mean(flux_data_to_avrg[:,:], axis=0)#time_average(fluxdata, interval=time_avrg)
        fluxes_err  = std(flux_data_to_avrg[:,:], axis=0)


    if save == True:
        fn = outdir(case) + '/fluxes.h5'
        hf = h5py.File(fn, 'w')
        hf.create_dataset('nspec', data=nspec(case))
        hf.create_dataset('fprim', data=nprim(case))
        hf.create_dataset('tprim', data=tprim(case))
        hf.create_dataset('ref_values', data=ref_values(case))
        hf.create_dataset('fluxes', data=fluxdata)
        if trange:
            hf.create_dataset('trange', data=trange)
            hf.create_dataset('fluxes_avrg', data=fluxes_avrg)
            hf.create_dataset('fluxes_err', data=fluxes_err)
        hf.close()
        print("File saved: ", fn)
    if plot:
        if si:
            xlab, ylab = '$t$ [$s$]',\
                         ['$\\Gamma_s$ [m$^{-2}$s$^{-1}$] ',\
                          '$\\Pi_s$ [kg s$^{-1}$]', '$Q_s$ [W m$^{-2}$]']
        else:
            xlab, ylab = '$t$ [$\\Omega_{r}^{-1}$]',\
                         ['$\\Gamma_s/\\Gamma_{gBs}$', '$\\Pi_s/\\Pi_{gBs}$', '$Q_s/Q_{gBs}$']
            
        for j in arange(0,3):
            ax = pl2d(xrange=None, yrange=None, xlabel=xlab, ylabel=ylab[j],\
                      fig_size=(8.5, 7.5))
            
            #for k in arange(nspecies):
            xran = [min(fluxdata[:,0]), max(fluxdata[:,0])]
            yran = [min(fluxdata[:,1+j*nspecies+spec])*1.05, max(fluxdata[:,1+j*nspecies+spec])*1.05]

            ax.plot(fluxdata[:,0], fluxdata[:,1+j*nspecies+spec],\
                    '-', color='black', linewidth=3, label='species = ' + str(spec))
            ax.set_xlim(xran)
            ax.set_ylim(yran)
            ax.ticklabel_format(style='sci', scilimits=(0,0))
            
            if trange:
                ax.axhline(y=fluxes_avrg[j*3+spec],\
                           linestyle='--', linewidth=2, color='blue')
                
                ax.barh(y=yran[0], width=(trange[1]-trange[0]),\
                        height=yran[1]-yran[0], \
                        left=trange[0], align='edge', facecolor='yellow', alpha=0.3)
                
        ax.legend(loc='best',labelspacing=0.0, prop={'size':26})
            
        show()

    if trange : return(fluxes_avrg, fluxes_err)
    else: return(fluxdata)


def time_average(data, time_column=0, interval=None):
    # Assuming the data is stored in an array such that
    # time is stored at data[:,time_average], this function averages
    # in in interval [t0,tf] the data store in data[:,j]
    # with j different from time_average.
    #
    # interval := [t0, tf]
    #
    data_avrg=empty(size(data[0,:]), dtype='float')
    data_avrg[:]=0
    num_val=0

    
    for i in arange(0, size(data[:, time_column])):
        if data[i, time_column] > interval[0] and data[i, time_column] < interval[1]:
            data_avrg[:]=data_avrg[:]+data[i,:]
            num_val = num_val+1

    return data_avrg/num_val
        



def ql_fluxes(case=None, view=False, last=True, formula='per'):
    i        = case
    tfull    = time(i)[0][1:]
    omegafile= outfile(i,'omega')
    fluxfile = outfile(i,'fluxes')
    n_time   = time(i)[1] # Including t=0
    n_ky     = ky(i)[1]
    n_kx     = kx(i)[1]
    ns       = nspec(i)
    density  = dens(i)[0]/dens(i)[0][0]
    # Full vectores to calculate the flux 
    flux    = empty((ns*3,size(tfull)), dtype='float')
    flux_ql = empty((ns*3,size(tfull)), dtype='float')
    flux_ql_last = empty(ns*3, dtype='float')
        
    fluxcols = ns*3+1

    if n_ky != 1 or n_kx != 1:
        print("Error: to use this function n_ky or n_kx cannot be > 1."); return

    omegadata = loadtxt(omegafile, dtype='float').reshape(n_time-1, n_kx, n_ky, 7)
    omega     = omegadata[:,0,0,3] 
    gamma     = omegadata[:,0,0,4]
    # First line of fluxdata is at time=0.0, omegadata begins after first time step
    print(fluxfile)
    fluxdata  = loadtxt(fluxfile, dtype='float')[1:,:]
    phi2data  = phi2_vs_kxky(i)[1:,0,0]
    kxrhoi    = kx(i)[0][0]
    kyrhoi    = ky(i)[0][0]
    
    # For plotting the particle fluxes
    if view:
        ax1 = pl2d(xrange=None, yrange=None, xlabel="$t v_{th,i}/a$",\
                   ylabel="$\\Gamma/\\Gamma_{gB}$", fig_size=(8.5, 7.5),\
                   title="$(k_x, k_y)\\rho_i=("+str(kx(i)[0][0])+", " + str(ky(i)[0][0])+")$",\
                   log=False)
        line = ['-', '--', ':']
        color = ['blue', 'crimson', 'green']

    for i in arange(0, 3*ns):
        flux[i,:]    = fluxdata[:,i+1]
        dens_i       = i-int(i/ns)*ns
        if formula == 'mik':
            flux_ql[i,:] = flux[i,:]*gamma/phi2data/density[dens_i]/kyrhoi**2.0
        elif formula == 'per':
            flux_ql[i,:] = flux[i,:]*(gamma**2.0+omega**2.0)/phi2data/density[dens_i]/kyrhoi**2.0/gamma
        # We remove not finite values
        
        notnan_flux      = isfinite(flux_ql[i,:])
        notnan_omega     = isfinite(omega[:])

        # We prepare all quantities to build the QL fluxes
        # up to the last value where these are finite
        # It can occur that this happens before than the instant
        # where omega finds the same problem.
        t_f        = tfull[notnan_flux]
        flux_ql_f  = flux_ql[i, notnan_flux]
        flux_f     = flux[i, notnan_flux]
        omega_ql   = omega[notnan_flux]
        gamma_ql   = gamma[notnan_flux]
        phi2data_f = phi2data[notnan_flux]

        omega_full = omega[notnan_omega]
        gamma_full = gamma[notnan_omega]
        omega_full_last = omega_full[size(omega_full)-1]
        gamma_full_last = gamma_full[size(gamma_full)-1]
        
        t_last_pos = size(t_f)-2
        
        t_last          = t_f[t_last_pos]
        flux_ql_last[i] = flux_ql_f[t_last_pos]
        omega_ql_last   = omega_ql[t_last_pos]
        gamma_ql_last   = gamma_ql[t_last_pos]
        phi2data_last   = phi2data[t_last_pos]
        
        if view and i < ns:
            ax1.plot(t_f, flux_ql_f, line[i], color=color[i], linewidth=4, label='species ='+str(i))

    if view:
        ax1.autoscale()
        ax1.legend(loc='best',labelspacing=0.0, prop={'size':24})
        show()

    data = list([kxrhoi, kyrhoi]) + list([omega_ql_last, gamma_ql_last, omega_full_last, gamma_full_last, phi2data_last]) +\
           list(flux_ql_last)
    # We return the number of species to know how many columns are expected
    return ns, array(data)
    
def phi(case, yrange=None, xrange=None):
    datafile    = outfile(case, quant='final_fields')
    
    vzed, nzed, iz0  = zed(case)[0], zed(case)[1], zed(case)[2]
    n_ky        = ky(case)[1]
    n_kx        = kx(case)[1]
    data        = loadtxt(datafile, dtype='float').reshape(nzed, n_kx, n_ky, 9)

    if n_kx == 1 and n_ky == 1:
        phi_re  = data[:, 0,0,4]
        phi_im  = data[:, 0,0,5]

        xlims = [vzed.min()*pi, vzed.max()*pi]
        ylims = [min(phi_im.min(), phi_re.min()), max(phi_im.max(), abs(phi_re).max())]
        
        pl2y(vzed*pi, phi_re, phi_im, xlabel='$\\zeta$',\
             ylabel='$\\tilde{\\varphi}(\\zeta)$',\
             title="$(k_x, k_y)=("+str(kx(case)[0][0])+", " + str(ky(case)[0][0])+")$",\
             xrange=xlims, yrange=ylims,\
             key1='$\\Re(\\tilde{\\varphi})(\\zeta)$',\
             key2='$\\Im(\\tilde{\\varphi})(\\zeta)$',\
             ax=None)        
        show()


def get_final_fields_txt(case):
    datafile=outfile(case, quant='final_fields')
    vzed, nzed, iz0  = zed(case)[0], zed(case)[1], zed(case)[2]
    n_ky        = ky(case)[1]
    n_kx        = kx(case)[1]
    data        = loadtxt(datafile, dtype='float').reshape(nzed, n_kx, n_ky, 11)
    phi2     = data[:,:,:,4]**2.0+data[:,:,:,2]**2.0
    phi2max  = phi2.max()
    phi2norm = phi2/phi2max
    
    return(phi2norm)

def phi_last(case, ikx=0, iky=0):
    phi                      = phi_vs_t(case)[0]
    nt, nzed, nkx, nky, nphi = shape(phi)
    vzed, nzed, iz0          = zed(case)[0], zed(case)[1], zed(case)[2]
    
    phi_last = phi[nt-1,:,ikx,iky,:]
    phi_re   = phi_last[:,0]
    phi_im   = phi_last[:,1]
    
    xlims = [vzed.min(), vzed.max()]
    ylims = [min(phi_im.min(), phi_re.min()), max(phi_im.max(), abs(phi_re).max())]
        
    pl2y(vzed, phi_last[:,0], phi_last[:,1], xlabel='$\\zeta$',\
         ylabel='$\\tilde{\\varphi}(\\zeta)$',\
         title="$(k_x, k_y)=("+str(kx(case)[ikx][iky])+", " + str(ky(case)[ikx][iky])+")$",\
         xrange=xlims, yrange=ylims,\
         key1='$\\Re(\\tilde{\\varphi})(\\zeta)$',\
         key2='$\\Im(\\tilde{\\varphi})(\\zeta)$',\
         ax=None)
    show()
    

def vmecgeo(case):
    # Function to read and represent geometric quantities

    zeta  = zed(case)[0]
    quant = geo(case)


    # Labels
    bmag     = '$B/B_{\mathrm{ref}}$'
    gradpar  = '$L_{\\mathrm{ref}}\\nabla_{\|} z$'
    gbdrift  = '$2B_{\\mathrm{ref}}L_{\mathrm{ref}}\\mathbf{B}\\times\\nabla B\\cdot\\nabla y/B^3$'
    gbdrift0 = '$\\hat{s}2B_{\mathrm{ref}}L_{\mathrm{ref}}\\mathbf{B}\\times\\nabla B\\cdot\\nabla x/B^3$'
    cdrift   = '$2B_{\mathrm{ref}}L_{\mathrm{ref}}\\mathbf{B}\\times\mathbf{\kappa}\\cdot\\nabla y/B^2$'
    cdrift0  = '$\\hat{s}2B_{\mathrm{ref}}L_{\mathrm{ref}}\\mathbf{B}\\times\\kappa\\cdot\\nabla x/B^2$'

    # Figures and axes
    
    fig = plt.figure(figsize=(18, 12))
    fig.subplots_adjust(left=0.1, wspace=0.35)
    
    ax1 = fig.add_subplot(231)
    ax2 = fig.add_subplot(232)
    ax3 = fig.add_subplot(233)
    ax4 = fig.add_subplot(234)
    ax5 = fig.add_subplot(235)
    ax6 = fig.add_subplot(236)

    plxy(zeta, quant[0], xlabel='$\\zeta$', ylabel=bmag,      title='\\texttt{bmag}'   , ax=ax1)
    plxy(zeta, quant[1], xlabel='$\\zeta$', ylabel=gradpar,   title='\\texttt{gradpar}', ax=ax2)
    plxy(zeta, quant[2], xlabel='$\\zeta$', ylabel=gbdrift,   title='\\texttt{gbdrift}', ax=ax3)
    plxy(zeta, quant[3], xlabel='$\\zeta$', ylabel=gbdrift0,  title='\\texttt{gbdrift0}',ax=ax4)
    plxy(zeta, quant[4], xlabel='$\\zeta$', ylabel=cdrift,    title='\\texttt{cdrift}',  ax=ax5)
    plxy(zeta, quant[5], xlabel='$\\zeta$', ylabel=cdrift0,   title='\\texttt{cdrift0}', ax=ax6)    

    fig2 = plt.figure(figsize=(18, 6))
    fig2.subplots_adjust(left=0.1, wspace=0.35)
    
    ax7 = fig2.add_subplot(131)
    ax8 = fig2.add_subplot(132)
    ax9 = fig2.add_subplot(133)

    gs2   = '$|\\nabla y|^2$'
    gds21 = '$\hat{s}^2\\nabla x\\cdot\\nabla y$'
    gds22 = '$\hat{s}^2|\\nabla x|^2$' 
    
    plxy(zeta, quant[6], xlabel='$\\zeta$', ylabel=gs2,      title='\\texttt{gds2}'   , ax=ax7)
    plxy(zeta, quant[7], xlabel='$\\zeta$', ylabel=gds21,    title='\\texttt{gds21}',   ax=ax8)
    plxy(zeta, quant[8], xlabel='$\\zeta$', ylabel=gds22,    title='\\texttt{gds22}',   ax=ax9)    


    show()

    
#    bmag, gradpar, gbdrift, gbdrift0, cvdrift, cvdrift0, gds2, gds21, gds22
    

def omega_max(case):
    return omega(case, last=True)[:,:,3].max(), omega(case, last=True)[:,:,4].max()


def ktkn(hpc=None, equil=None, process=False):

    hpc ='marconis'
    equil = 'w7xr003s'
    nLn = 21
    nLT = 21
    nky = 20
    run0=101
    irun=run0
    outfile=datadir() + 'marconis_w7xr003s_0101_0441.dat'

    omega_r     = empty((nLn, nLT, nky), dtype='float')
    omega_i     = empty((nLn, nLT, nky), dtype='float')
    omega_r_max = empty((nLn, nLT),      dtype='float')
    omega_i_max = empty((nLn, nLT),      dtype='float')    
    mLn         = empty((nLn, nLT),      dtype='float')
    mLT         = empty((nLn, nLT),      dtype='float')

    if process == True:
        # The runs are analyzed and the file compiling all the data
        # is written
        f = open(outfile, 'w')
        f.write('# (0) run   (1) i_Ln  (2) i_LT   (3) nprim   (4) tprim   '+\
                '(5) Max(Re(omega)) (6) Max(Im(omega))' + '\n')
        
        print("Writing file :", outfile + '\n')
        
        for i in arange (0, nLn):
            print('# (0) run   (1) i_Ln  (2) i_LT   (3) nprim   (4) tprim   '+\
                  '(5) Max(Re(omega)) (6) Max(Im(omega))')
            
            for j in arange (0, nLT):
                
                if hpc != None:
                    case  = hpc + '/' + equil + '_' + str(format8(irun))
                else:
                    case  = equil + '_' + str(format8(irun)) 
                    
                nprim = read_stella_float(case, 'fprim')[0]
                tprim = read_stella_float(case, 'tprim')[0]

                mLn[i, j]         = nprim
                mLT[i, j]         = tprim
                
                omega_r[i, j, :]  = omega(case, last=True)[0, :, 3]
                omega_i[i, j, :]  = omega(case, last=True)[0, :, 4]
                omega_r_max[i, j] = omega(case, last=True)[0, :, 3].max()
                omega_i_max[i, j] = omega(case, last=True)[0, :, 4].max()
                
                print(case, format8(i), format8(j), format2(nprim), format2(tprim),\
                      format2(omega_r[i, j, :].max()), format2(omega_i[i, j, :].max()))
                
                f.write(case+'\t'+str(format8(i))+'\t'+str(format8(j))+'\t'+\
                        str(format2(nprim))+'\t'+str(format2(tprim))+'\t'+\
                        str(format2(omega_r[i, j, :].max()))+'\t'+\
                        str(format2(omega_i[i, j, :].max()))+'\n')
                          
                irun = irun + 1

        f.close()

    if process == False:
        # Just read the existing file with all the data from the set of runs.
        print("Reading existing file : ", outfile)

        data = loadtxt(outfile, dtype='float', usecols=[3,4,5,6])
        mLn  = data[:,0].reshape((nLn, nLT))
        mLT  = data[:,1].reshape((nLn, nLT))
        omega_r_m = data[:,2].reshape((nLn, nLT))
        omega_i_m = data[:,3].reshape((nLn, nLT))      

        fig, ax = plt.subplots(figsize=(11.0, 8.4))
        col = ax.pcolormesh(mLn[:,0], mLT[0,:], omega_r_m.T, vmin=0, vmax=1.8)
#        CS1 = ax.contour((mLn[:,0], mLT[0,:], omega_r_m.T, levels=[0.1, 0.3, 0.5, 0.7, 0.9, 1.0], colors='white')
#        ax.clabel(CS1, inline=1, fontsize=14,  fmt='%1.1f', linewidths = 1.5)
        ax.set_xlabel('$a/L_{n_i}$')
        ax.set_ylabel('$a/L_{T_i}$')
        ax.xaxis.set_minor_locator(AutoMinorLocator(10))
        ax.yaxis.set_minor_locator(AutoMinorLocator(10))
        cbar=plt.colorbar(col)
        cbar.set_label('$\\Im({\\omega})^{\\textrm{max}}$')
        cbar.set_label('$a\\Im({\\omega})^{max}/v_{\\mathrm{th,i}}$')
        ax.set_aspect(1.0)
        title = '$k_{x}=0$; $0\\le k_{y}\le 10$'
        ax.set_title(title, y=1.02)
        ax.grid(color='grey', linestyle='-', linewidth=0.3)
                
                #        print mLn[:,0]
                #        print mLT[0,:]
                #        print omega_i_m
                #        mLT  = data[0,:].reshape(nLn, nLT)
                
                #    cmap(xdata=mLn[:,0], ydata=mLT[0,:], zdata=omega_r_m,\
                #         xlabel='$-a/L_{n_i}$', ylabel='$-a/L_{T_i}$', zlabel='$\\Im(\\omega)/v_{t}$')
                
        show()

def rephi_vs_kxky(case):
    # electrostatic real potential averaged as function of (ky,kx,t)
    rephi_vs_kxky_stella = read_stella_float(case, 'phi_real')
    return rephi_vs_kxky_stella

def imphi_vs_kxky(case):
    # electrostatic real potential averaged as function of (ky,kx,t)
    imphi_vs_kxky_stella = read_stella_float(case, 'phi_imag')
    return imphi_vs_kxky_stella

def ZF(case, plot=False, verbose=True, fit=True):
    #It computes the ZF response using <Re(phi)>:
    t = time(case)[0]
    re_phi_avg = rephi_vs_kxky(case)[:,0,0]
    im_phi_avg = imphi_vs_kxky(case)[:,0,0]
    re_phi_avg_norm = re_phi_avg/max([max(re_phi_avg), max(im_phi_avg)])
    im_phi_avg_norm = im_phi_avg/max([max(re_phi_avg), max(im_phi_avg)])

    if fit:
        popt, _ = curve_fit(objective, t, re_phi_avg_norm)

        if verbose:
            print('A_ZF     [addim.]   = ', str(format2(popt[0])))
            print('Omega_ZF [v_ref/a]  = ', str(format2(popt[1])))
            print('gamma_ZF [v_ref/a]  = ', str(format2(popt[2])))
            print('R_ZF     [adim.]    = ', str(format2(popt[3])))
            print('c_1      [adim.]    = ', str(format2(popt[4])))
            print('c_2      [adim.]    = ', str(format2(popt[5])))
            print('c_3      [a.u.]     = ', str(format2(popt[6])))


    # Writing out data to txt file ==============================================
    outfile = runsdir() + '/' + case + '/ZF_stella.dat'
    f = io.open(outfile, 'w')

    f.write('# Fitting function for <Re[phi(t)]/Re[phi(t=0)]>:\n')
    if fit:
        f.write('# f(t)=A_ZF * cos(Omega_ZF*t) * exp(-gamma_ZF*t) + R_ZF + c_1/(c_2+t**c_3)\n')
        f.write('# (0)A_ZF[addim.]\t (2)Omega_ZF[v_ref/a]\t (3)gamma_ZF[v_ref/a]\t (4)R_ZF[addim.]\t (4)c_1[addim.]\t (5)c_2[addim.]\t (6)c_3[a.u.]\n')
        f.write('#'+str(format2(popt[0]))+'\t'+str(format2(popt[1]))+'\t'+str(format2(popt[2]))+'\t'+\
                str(format2(popt[3]))+'\t'+str(format2(popt[4]))+'\t'+str(format2(popt[5]))+'\t'+\
                str(format2(popt[6]))+'\n')
    f.write('# (0) t [a/vref]   (1) <Re(phi)>/<Re(phi)>_max   (2) <Im(phi)>/<Im(phi)>_max' + '\n')
    
    for i in arange(0, size(t)):
        f.write(str(format2(t[i]))+'\t'+str(format2(re_phi_avg_norm[i]))+'\t'+str(format2(im_phi_avg_norm[i]))+'\n')
    f.close()
    print("File saved: ", outfile)
    
    # Writing out data to h5 file ==============================================
    outfile = runsdir() + '/' + case + '/ZF_stella.h5'
    f = h5py.File(outfile, 'w')
    f.create_dataset('time', data=t)
    f.create_dataset('re_phi_avr_norm', data=re_phi_avg_norm)
    f.create_dataset('im_phi_avg_norm', data=im_phi_avg_norm)
    if fit:
        f.create_dataset('Fit', data='A_ZF * cos(Omega_ZF*t) * exp(-gamma_ZF*t) + R_ZF + c_1/(c_2+t**c_3)')    
        f.create_dataset('A_ZF', data=popt[0])
        f.create_dataset('Omega_ZF', data=popt[1])
        f.create_dataset('gamma_ZF', data=popt[2])
        f.create_dataset('R_ZF',     data=popt[3])
        f.create_dataset('c_1',      data=popt[4])
        f.create_dataset('c_2',      data=popt[5])
        f.create_dataset('c_3',      data=popt[6])    
        f.close()
    print("File saved: ", outfile)

    if plot:
        ax= pl2d(xlabel='$t v_{th}/\\rho_s$', xrange=[0,max(t)],\
                 yrange=[-1,1],\
                 ylabel='$\\langle\\Re(\\varphi)\\rangle/\\max(\\langle\\Re(\\varphi)\\rangle)$')
        ax.plot(t,re_phi_avg_norm,linestyle='-', color='k', linewidth=2)
        ax.plot(t,im_phi_avg_norm,linestyle='-', color='r', linewidth=1)        

    show()


    
    return t, re_phi_avg_norm

# Objective function

def objective(t, A_ZF, Omega_ZF, gamma_ZF, R_ZF, c_1, c_2, c_3):

    return A_ZF * cos(Omega_ZF*t) * exp(-gamma_ZF*t) + R_ZF + c_1/(c_2+t**c_3)
    
