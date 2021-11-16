
#================================================================
# Plot |phi(kx,ky)^2| averaged over the time range with steady state fluxes
#================================================================

import os, sys
import numpy as np
from scipy.interpolate import interp2d
# from .plotbox import plotbox_2d
# import matplotlib.gridspec as gridspec
# import matplotlib.pyplot as plt
from stellapy.utils.decorators import verbose_wrapper
# from stellapy.utils import ensure_dir, get_PathNameOfFolder, get_NetcdfFilesOrReducedFilesInside
# from stellapy.decorators import verbose_wrapper, exit_program

@verbose_wrapper
def plot_potentialVsKxVsKy(
        # Get the data
        folder, \
        # Axes ranges
        t_start=400, \
        t_end=None, \
        t_range=None, \
        c_range=None, \
        log=True, \
        # Smooth the image
        interpolate=4, \
        # Toggles
        plot_fluxes=False, \
        plot_omega=False, \
        plot_gamma=False, \
        compare_radii=False, \
        time_evolution=True):
    ''' Plot |phi(kx,ky)^2| averaged over the time range with steady state fluxes. ''' # t_start=300, t_end=500

    from stellapy.decorators.verbose_wrapper import indent    

    # Choose the plotting option
    if compare_radii: log=False

    # Initiate the folders
    folders = [folder]

    # Initiate figures
    if not compare_radii and not time_evolution and not plot_fluxes:
        show_colorbar, fig, ax = initiate_figures(folder, compare_radii, time_evolution, plot_fluxes)
    if time_evolution:
        show_colorbar, fig, ax1,ax2,ax3 = initiate_figures(folder, compare_radii, time_evolution, plot_fluxes)
    if plot_fluxes and not compare_radii: 
        show_colorbar, fig, ax1,ax2,ax3,ax11,ax12,ax21,ax22,ax31,ax32 = initiate_figures(folder, compare_radii, time_evolution, plot_fluxes)
    if compare_radii:
        show_colorbar, fig, ax1,ax2,ax3,ax4 = initiate_figures(folder, compare_radii, time_evolution, plot_fluxes)       
        folders = ["nl_w7xr169+252_0001","nl_w7xr169+252_0002","nl_w7xr169+252_0003_y015","nl_w7xr169+252_0004_y015"]   
        folders = ["nl_w7xr169+252_0001","nl_w7xr169+252_0002","nl_w7xr169+252_0003","nl_w7xr169+252_0004"]         
        folders = ["nl_w7xr169+252_0001_y015","nl_w7xr169+252_0002_y015","nl_w7xr169+252_0003_y015","nl_w7xr169+252_0004_y015"]
        folders = ["nl_w7xr348+252_0001", "nl_w7xr348+252_0002", "nl_w7xr348+252_0003_broken", "nl_w7xr348+252_0004"] 
        folders = ["nl_w7xr348+252_0001_y015", "nl_w7xr348+252_0002_y015", "nl_w7xr348+252_0003_y015", "nl_w7xr348+252_0004_y015"] 
        folders = ["nl_W7xr348+s17_0001_y015","nl_W7xr348+s17_0002_y015","nl_W7xr348+s17_0003_y015","nl_W7xr348+s17_0004_y015"] 
        folders = ["nl_W7xr348+s17_0001","nl_W7xr348+s17_0002","nl_W7xr348+s17_0003","nl_W7xr348+s17_0004"] 

    # Iterate over casefolders
    for folder in folders:

        # Get all out.h5 files OR out.nc files inside <folder>
        netcdf_files = get_NetcdfFilesOrReducedFilesInside(folder)

        # Read the netcdf files to get potential['time', 'phi2_vs_kxky', 'phi2'] or potential[ky]['time', 'phi2_vs_kxky', 'phi2']
        potential = read.phi(folder)
        keys_shot, keys_rho = read.sorted_shot_and_rho_keys(potential, return_type='keys')

        # Read the data
        shot     = keys_shot[0]
        rho      = keys_rho[shot][0]
        ky       = list(potential[shot][rho].keys())[0]
        vec_time     = potential[shot][rho][ky]['time']
        phi2         = potential[shot][rho][ky]['phi2']
        phi_log      = np.log10(np.abs(phi2))
        phi2_vs_kxky = potential[shot][rho][ky]['phi2_vs_kxky']
        phi2_vs_t    = potential[shot][rho][ky]['phi2_vs_t']        # phi(t,kx,ky,ri)
        dim_kx       = potential[shot][rho][ky]['dim_kx']
        vec_kx       = potential[shot][rho][ky]['vec_kx'] 
        dim_ky       = potential[shot][rho][ky]['dim_ky']
        vec_ky       = potential[shot][rho][ky]['vec_ky']  
        if plot_fluxes and 'time_fluxes' in potential[shot][rho][ky]:
            vec_time     = potential[shot][rho][ky]['time_fluxes'] 
            qflx_vs_kxky = np.abs(potential[shot][rho][ky]['qflx_kxky'])
            pflx_vs_kxky = np.abs(potential[shot][rho][ky]['pflx_kxky'])
            vflx_vs_kxky = np.abs(potential[shot][rho][ky]['vflx_kxky'])
            qflx         = np.abs(potential[shot][rho][ky]['qflx'])
            pflx         = np.abs(potential[shot][rho][ky]['pflx'])
            vflx         = np.abs(potential[shot][rho][ky]['vflx'])
        if plot_fluxes and not 'time_fluxes' in potential[shot][rho][ky]:
            print("EXIT: there are no fluxes found in the files.")
            return

        # Determine the time range to average over
        t_ranges = determine_trange(folder, compare_radii, time_evolution, plot_fluxes, vec_time, phi_log, t_range, t_start, t_end)

        #===================
        # Plot the data
        #===================

        # Add a plot of fluxes versus time
        if plot_fluxes and not time_evolution and not compare_radii:
            if '348' in folder: color='navy'
            if '169' in folder: color='crimson'
            ax11.plot(vec_time, qflx, lw=3., color=color, label="") 
            ax21.plot(vec_time, pflx, lw=3., color=color, label="") 
            ax31.plot(vec_time, vflx, lw=3., color=color, label="") 
        
        counter=-1
        for t_range in t_ranges:
            counter += 1

            # Initiate the axes
            vec_kx  = potential[shot][rho][ky]['vec_kx'] 
            vec_ky  = potential[shot][rho][ky]['vec_ky']  
            if not compare_radii: ax4 = None
            show_colorbar, ax, color, c_range, log = initiate_axes(folder, compare_radii, time_evolution, plot_fluxes, counter, ax1, ax2, ax3, ax4)
            if time_evolution: ax.set_title("$t\, v_{\mathrm{th}}/a=["+str(int(t_range[0]))+","+str(int(t_range[1]))+"]$", color=color)

            # Set the y-data    
            if plot_fluxes and compare_radii:  
                y_vs_kxky = qflx_vs_kxky; 
            elif plot_fluxes:       
                if counter==0:  y_vs_kxky = qflx_vs_kxky; 
                if counter==1:  y_vs_kxky = pflx_vs_kxky; 
                if counter==2:  y_vs_kxky = vflx_vs_kxky;  
            elif plot_fluxes and time_evolution:       
                if counter==0:  y_vs_kxky = qflx_vs_kxky; 
                if counter==1:  y_vs_kxky = qflx_vs_kxky; 
                if counter==2:  y_vs_kxky = qflx_vs_kxky;  
            elif not plot_fluxes: y_vs_kxky = phi2_vs_kxky;             
            elif plot_omega:      y_vs_kxky = np.log10(phi2_vs_t[:,:,:,:,0]);    
            elif plot_gamma:      y_vs_kxky = np.log10(phi2_vs_t[:,:,:,:,1]); 
            
            if (plot_fluxes and t_range[0] >= vec_time[0] and t_range[1] <= vec_time[-1]) or not plot_fluxes: 

                # Calculate y_average over (kx,ky)
                vec_kx, vec_ky, y_avrg = average_y_over_kxky(y_vs_kxky, vec_time, t_range, interpolate, vec_kx, vec_ky)

                # Plot y_avrg(kx,ky)
                ax = cmap(xdata=vec_kx, ydata=vec_ky, zdata=y_avrg, xlabel='$k_x$', ylabel='$k_y$', zlabel='', title='', ax=ax, fig=fig,\
                          contour=1, cmap='jet', crange=c_range, log=log, compare_radii=compare_radii, show_colorbar=show_colorbar)

                if plot_fluxes and not time_evolution and not compare_radii:
                    if counter==0: ax = ax12; color='crimson'; c_range = None;  show_colorbar=True
                    if counter==1: ax = ax22; color='navy';    c_range = None;  show_colorbar=True
                    if counter==2: ax = ax32; color='green';   c_range = None;  show_colorbar=True 

                    ax = cmap(xdata=vec_kx, ydata=vec_ky, zdata=y_avrg, xlabel='$k_x$', ylabel='$k_y$', zlabel='', title='', ax=ax, fig=fig,\
                              contour=1, cmap='jet', crange=c_range, log=log, compare_radii=compare_radii, show_colorbar=show_colorbar)


        if plot_fluxes and not time_evolution and not compare_radii:
            for ax in [ax11, ax21, ax31]:
                ax.autoscale()

    plt.show()
    return


#================================================================
# Make a surface plot
#================================================================

# Plot a surface plot z(x,y) with a colormap
def cmap(xdata=None, ydata=None, zdata=None, xlabel='x', ylabel='y', zlabel='z', title=None,\
         ctics=False, num=None, contour=1, cmap='jet', vfield=False, vxdata=None, fig=None, show_colorbar=True,\
         vydata=None, cont2=None, epsname=None, crange=None, fig_size=(8.5, 7.5),ax=None, log=False, compare_radii=False):

    # Load the modules
    from matplotlib.colors import LogNorm
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    # For z(x,y) the columns refer to x and the rows to y
    zdata = np.transpose(zdata)

    # Value kx=0 is plotted from kx=0 to the next kx, correct this by shifting xdata half a tile left
    xdata = xdata-(xdata[-1]-xdata[-2])/2

    # Set the labels    
    ax.axis([xdata[0], xdata[np.size(xdata)-1], ydata[0], ydata[np.size(ydata)-1]])

    # Show nan values in black
    cmap = plt.get_cmap(cmap)
    cmap.set_bad(color='black')

    # Plot the surface plotz(x,y)
    if log == True:   norm = LogNorm(vmin=np.nanmin(zdata), vmax=np.nanmax(zdata))
    if log == False:  norm = None
    if not crange:    img  = ax.pcolormesh(xdata, ydata, abs(zdata), cmap=cmap, norm=norm)   
    if crange:        img  = ax.pcolormesh(xdata, ydata, zdata, cmap=cmap, vmin=crange[0], vmax=crange[1], norm=norm)
    
    # Add colorbar to specific sublpot
    if show_colorbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.1)
        cbar = plt.colorbar(img, cax=cax) 

#    if show_colorbar and log == False:  
#        #cbar.set_label(zlabel)
#        cbar.formatter.set_powerlimits((0,0))
#        cbar.formatter.set_scientific(True)
#        cbar.update_ticks()

    ax.set_ylim(ymin=0, ymax=2)
    fig.tight_layout(pad=0.15, h_pad=0, w_pad=0)
    return ax

#================================================================
# Initiate axes
#================================================================

def average_y_over_kxky(y_vs_kxky, vec_time, t_range, interpolate, vec_kx, vec_ky):
    from stellapy.decorators.verbose_wrapper import indent  
    # phi2_vs_kx_ky has dimensions (ntime, nakx, naky)
    # Here, we remove all values of phi2_vs_kx_ky of times out of the selected interval.
    y_to_avrg = y_vs_kxky[(vec_time > t_range[0]) & (vec_time < t_range[1]),:,:]
    
    # Each tile is the average of |phi(kx,ky)^2| over the time range
    y_avrg = np.nanmean(y_to_avrg, axis=0)

    # Remove broken tiles
    print(indent, 'Warning: all tiles with phi2<E-30 or phi2>100 are removed and shown in black in phi'+str(np.shape(y_avrg))+'.')
    print(indent, '   Tiles lower than E-30 are at:', np.argwhere(y_avrg[:,:]<1.E-30))
    print(indent, '   Their values are:', y_avrg[y_avrg[:,:]<1.E-30])
    print(indent, '   Tiles bigger than 100 are at:  ', np.argwhere(y_avrg[:,:]>100))
    print(indent, '   Their values are:', y_avrg[y_avrg[:,:]>100])
    filter_low  = y_avrg[:,:]<1.E-30 # kx,ky = 0,0 is 2.58485853e-32
    filter_high = y_avrg[:,:]>100 # kx,ky = left right,0 is 28
    
 
    # Interpolate the data
    if interpolate: 
        y_avrg[filter_low] = np.nan
        y_avrg[filter_high] = np.nan
        for index in np.argwhere(np.isnan(y_avrg)): # At ky,kx=0,0 the values are nan's
            if ( index[0]<=27 or index[0]==len(y_avrg[:,1])-1 ) and (not index[0]==0): 
                y_avrg[index[0], index[1]] = y_avrg[index[0]-1, index[1]]  
            if ( index[0]>=28 or index[0]==0 ) and (not index[0]==len(y_avrg[:,1])-1):           
                y_avrg[index[0], index[1]] = y_avrg[index[0]+1, index[1]]
        function = interp2d(vec_ky, vec_kx, y_avrg, kind='linear')
        xnew = np.linspace(vec_ky[0], vec_ky[-1], int(len(vec_ky))*interpolate)
        ynew = np.linspace(vec_kx[0], vec_kx[-1], int(len(vec_kx))*interpolate)
        y_avrg_interp = function(xnew,ynew)
        vec_ky_interp, vec_kx_interp = np.meshgrid(xnew, ynew)
        vec_ky_interp = vec_ky_interp[0,:]; vec_kx_interp = vec_kx_interp[:,0]

    # Remove broken tiles
    if interpolate: 
        y_avrg_interp[np.repeat(np.repeat(filter_low,  interpolate, axis=1), interpolate, axis=0)] = np.nan
        y_avrg_interp[np.repeat(np.repeat(filter_high, interpolate, axis=1), interpolate, axis=0)] = np.nan
    if not interpolate:
        y_avrg[filter_low] = np.nan
        y_avrg[filter_high] = np.nan

    if not interpolate:     return vec_kx, vec_ky, y_avrg
    if interpolate:         return vec_kx_interp, vec_ky_interp, y_avrg_interp

#================================================================
# Initiate axes
#================================================================

def initiate_axes(folder, compare_radii, time_evolution, plot_fluxes, counter, ax1, ax2, ax3, ax4):

    log=True
    if '348' in folder: color='navy'
    if '169' in folder: color='crimson'
    if '+s17' in folder: color='orange'

    if compare_radii:
        log = True
        c_range = [1.E-2,1E-1]
        c_range = [0,0.01]
        c_range = [1.E-3,1]
        c_range = [1.E-2,10]
        c_range = [1.E-3,1]
        if plot_fluxes: c_range = [0,0.01]
        if plot_fluxes and log==True: c_range = [1.E-3,1.E-1]
        if "001" in folder:    ax = ax1;   show_colorbar=False
        if "002" in folder:    ax = ax2;   show_colorbar=False
        if "003" in folder:    ax = ax3;   show_colorbar=False
        if "004" in folder:    ax = ax4;   show_colorbar=True
    if time_evolution and not compare_radii:
        if counter==0:  ax = ax1;    color='crimson'; c_range = [1.E-10,1E-1];  show_colorbar=True; log=True
        if counter==0:  ax = ax1;    color='crimson'; c_range = [1.E-6,1E-3];  show_colorbar=True; log=True
        if counter==1:  ax = ax2;    color='navy';    c_range = [1.E-2,10];   show_colorbar=False; log=True
        if counter==2:  ax = ax3;    color='green';   c_range = [1.E-2,10];   show_colorbar=True; log=True
#        if counter==1:  ax = ax2;    color='navy';    c_range = [1.E-2,1E-1];   show_colorbar=False; log=False
#        if counter==2:  ax = ax3;    color='green';   c_range = [1.E-2,1E-1];   show_colorbar=True; log=False
#        if counter==0:  ax = ax1;    color='crimson'; c_range = None;           show_colorbar=True
#        if counter==1:  ax = ax2;    color='navy';    c_range = None;           show_colorbar=True
#        if counter==2:  ax = ax3;    color='green';   c_range = None;           show_colorbar=True
    if plot_fluxes and not compare_radii:
        if counter==0:  ax = ax1;    color='crimson'; c_range = [1.E-3,2E-2];    show_colorbar=True
        if counter==1:  ax = ax2;    color='navy';    c_range = [8.E-13,1E-11];  show_colorbar=True
        if counter==2:  ax = ax3;    color='green';   c_range = [1.E-5,2.E-4];   show_colorbar=True 

                    
    return show_colorbar, ax, color, c_range, log

#================================================================
# Determine time range
#================================================================

def determine_trange(folder, compare_radii, time_evolution, plot_fluxes, vec_time, phi_log, t_range, t_start, t_end):
    from stellapy.decorators.verbose_wrapper import indent  
    # If the start of the time range is given, the end is the last time value or t_end  
    if plot_fluxes==True:                       t_range =[np.max([400.0,vec_time[0]]), vec_time[-1]] # Plot over the saturated phase
    if plot_fluxes==True and "003" in folder: t_range[0]=t_range[0]-100
    if plot_fluxes==True and "004" in folder: t_range[0]=t_range[0]-200
    if not t_range and t_start and not t_end:   t_range = [t_start, vec_time[-1]]
    if not t_range and t_start and t_end:       t_range = [t_start, np.min([vec_time[-1], t_end])]
    if t_range and t_range[1] < t_range[0]:
        exit_program('There is no data in the chosen time range', phi2_vs_kxky_averaged, sys._getframe().f_lineno)
    print(" The time axis range is: ["+str(vec_time[0])+","+str(vec_time[-1])+"]")
    if plot_fluxes==False: print("     |phi(kx,ky)^2| is averaged over the time range: ["+str(t_range[0])+","+str(t_range[-1])+"]")
    if plot_fluxes==True:  print("     |Q/Q_{gB}^2| is averaged over the time range: ["+str(t_range[0])+","+str(t_range[-1])+"]")
    # Choose the time ranges
    if time_evolution:      
        # Show phi(kx,ky) in the linear growth
        time_start1 = vec_time[list(phi_log).index( list(filter(lambda i: i > -5.0, phi_log))[0] )]
        time_stop1  = vec_time[list(phi_log).index( list(filter(lambda i: i > -1.0, phi_log))[0] )]
        # Show phi(kx,ky) on the stagnation
        try: 
            time_start2 = vec_time[list(phi_log).index( list(filter(lambda i: i > 1.0, phi_log))[0] )]
            time_stop2  = vec_time[list(vec_time).index( list(filter(lambda i: i > time_start2+10, vec_time))[0] )]
        except:
            time_start2 = 10
            time_stop2  = 30
        # Show phi(kx,ky) in steady-state
        time_start3 = 400
        time_stop3  = 600 
        # Construct the t_ranges
        t_ranges = [[time_start1, time_stop1], [time_start2, time_stop2], [time_start3, time_stop3]]
    if not time_evolution:  t_ranges=[t_range]
    #elif not plot_fluxes:  t_ranges=[t_range]
    #elif plot_fluxes:      t_ranges=[t_range, t_range, t_range]
    return t_ranges
#================================================================
# Set up the figure
#================================================================

def initiate_figures(folder, compare_radii, time_evolution, plot_fluxes):

    # Set the labels 
    xlabel = '$k_x$' #'$k_{x}\\rho_i$'
    ylabel = '$k_y$' #'$k_{y}\\rho_i$'
    title  = '' #'$|\\varphi(k_x, k_y)|^2$'
    show_colorbar=True

    # Standard figure
    if not compare_radii and not time_evolution and not plot_fluxes:
        fig, ax = plt.subplots(figsize=(8.5, 7.5))
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)
        show_colorbar=True
        plot = None
        return show_colorbar, fig, ax

    # Compare radii 
    if compare_radii:
        plt.rc('font', size=20)
        fig = plt.figure(figsize=(15, 4))
        gs1 = gridspec.GridSpec(1, 100)
        gs1.update(wspace=0.025,bottom=0.2)
        x_label = "$k_{x}\\rho_i$"
        y_label = "$k_{y}\\rho_i$"  
        ax1 = plotbox_2d(x_label=x_label, y_label=y_label, title="$\\rho=0.5$", fig_size=(8.5, 7.5), ax=plt.subplot(gs1[0,:9]) );
        ax2 = plotbox_2d(x_label=x_label, y_label="", title="$\\rho=0.6$", fig_size=(8.5, 7.5), ax=plt.subplot(gs1[0,10:25]) );
        ax3 = plotbox_2d(x_label=x_label, y_label="", title="$\\rho=0.7$", fig_size=(8.5, 7.5), ax=plt.subplot(gs1[0,26:52]) );
        ax4 = plotbox_2d(x_label=x_label, y_label="", title="$\\rho=0.8$", fig_size=(8.5, 7.5), ax=plt.subplot(gs1[0,53:]) );
        ax2.yaxis.set_major_formatter(plt.NullFormatter())
        ax3.yaxis.set_major_formatter(plt.NullFormatter())
        ax4.yaxis.set_major_formatter(plt.NullFormatter())
        ax2.yaxis.set_major_formatter(plt.NullFormatter())
        ax3.yaxis.set_major_formatter(plt.NullFormatter())
        show_colorbar=False
        if "004" in folder: ax3.xaxis.set_major_formatter(plt.NullFormatter()); show_colorbar=True
        return show_colorbar, fig, ax1, ax2, ax3, ax4
    
    # Plot three subfigures: duringt the ramp-up of phi, during the stagnation and during steady-state
    if time_evolution:
        fig = plt.figure(figsize=(18, 5))
        gs1 = gridspec.GridSpec(1, 100)
        gs1.update(wspace=0.025,bottom=0.2)
        x_label = "$k_{x}\\rho_i$"
        y_label = "$k_{y}\\rho_i$"  
        ax1 = plotbox_2d(x_label=x_label, y_label=y_label, fig_size=(8.5, 7.5), ax=plt.subplot(gs1[0,:30]) );
        ax2 = plotbox_2d(x_label=x_label, y_label="", fig_size=(8.5, 7.5), ax=plt.subplot(gs1[0,39:69]) );
        ax3 = plotbox_2d(x_label=x_label, y_label="", fig_size=(8.5, 7.5), ax=plt.subplot(gs1[0,70:]) );
#        ax1 = plotbox_2d(x_label=x_label, y_label=y_label, fig_size=(8.5, 7.5), ax=plt.subplot(gs1[0,:28]) );
#        ax2 = plotbox_2d(x_label=x_label, y_label="", fig_size=(8.5, 7.5), ax=plt.subplot(gs1[0,34:63]) );
#        ax3 = plotbox_2d(x_label=x_label, y_label="", fig_size=(8.5, 7.5), ax=plt.subplot(gs1[0,70:]) );
        ax2.yaxis.set_major_formatter(plt.NullFormatter())
        ax3.yaxis.set_major_formatter(plt.NullFormatter())
        return show_colorbar, fig, ax1, ax2, ax3

    # Plot_fluxes
    if plot_fluxes and not time_evolution and not compare_radii:
        fig = plt.figure(figsize=(18, 5))
        gs = gridspec.GridSpec(1, 100)
        gs.update(wspace=0.025,bottom=0.2)
        x_label = "$k_{x}\\rho_i$"
        y_label = "$k_{y}\\rho_i$"
        ax1 = plotbox_2d(x_label=x_label, y_label=y_label, title = "$Q/Q_{gB}$", fig_size=(8.5, 7.5), ax=plt.subplot(gs[0,:28]) );
        ax2 = plotbox_2d(x_label=x_label, y_label="", title = "$\\Gamma/\\Gamma_{gB}$", fig_size=(8.5, 7.5), ax=plt.subplot(gs[0,34:62]) );
        ax3 = plotbox_2d(x_label=x_label, y_label="", title = "$\\Pi/\\Pi_{gB}$", fig_size=(8.5, 7.5), ax=plt.subplot(gs[0,69:]) );
        ax2.yaxis.set_major_formatter(plt.NullFormatter())
        ax3.yaxis.set_major_formatter(plt.NullFormatter())
        fig1 = plt.figure(figsize=(12, 7))
        gs1 = gridspec.GridSpec(1, 2)
        gs1.update(wspace=0.025,bottom=0.2)
        ax11 = plotbox_2d(x_label="$t\, v_{\mathrm{th}}/a$", y_label="$Q/Q_{gB}$", fig_size=(8.5, 7.5), ax=plt.subplot(gs1[0]));
        ax12 = plotbox_2d(x_label="$k_{x}\\rho_i$", y_label="$k_{y}\\rho_i$", fig_size=(8.5, 7.5), ax=plt.subplot(gs1[1]));
        fig2 = plt.figure(figsize=(12, 6))
        gs2 = gridspec.GridSpec(1, 2)
        gs2.update(wspace=0.025,bottom=0.2)
        ax21 = plotbox_2d(x_label="$t\, v_{\mathrm{th}}/a$", y_label="$\\Gamma/\\Gamma_{gB}$", fig_size=(8.5, 7.5), ax=plt.subplot(gs2[0]));
        ax22 = plotbox_2d(x_label="$k_{x}\\rho_i$", y_label="$k_{y}\\rho_i$", fig_size=(8.5, 7.5), ax=plt.subplot(gs2[1]));
        fig3 = plt.figure(figsize=(15, 6))
        gs3 = gridspec.GridSpec(1, 2)
        gs3.update(wspace=0.025,bottom=0.2)
        ax31 = plotbox_2d(x_label="$t\, v_{\mathrm{th}}/a$", y_label="$\\Pi/\\Pi_{gB}$", fig_size=(8.5, 7.5), ax=plt.subplot(gs3[0]));
        ax32 = plotbox_2d(x_label="$k_{x}\\rho_i$", y_label="$k_{y}\\rho_i$", fig_size=(8.5, 7.5), ax=plt.subplot(gs3[1]));
        return show_colorbar, fig, ax1, ax2, ax3, ax11, ax12, ax21, ax22, ax31, ax32
