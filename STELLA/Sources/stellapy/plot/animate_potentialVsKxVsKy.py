
#================================================================
# Animate |phi(kx,ky)^2| in function of time
#================================================================

import numpy as np
from scipy.interpolate import interp2d
from stellapy.utils.decorators import verbose_wrapper
# from stellapy.utils import get_PathNameOfFolder, get_NetcdfFilesOrReducedFilesInside

@verbose_wrapper # stellapy.plot.kspectra_movie
def animate_potentialVsKxVsKy(\
        # Get the data
        folder, \
        # Axes ranges
        log=True, \
        t_range=None, \
        # Change number of frames and the smoothness
        time_step=0.5, \
        interpolate=4, \
        # Color of the plots
        cmap='jet', \
        # Toggles
        plot_phi=True, \
        plot_start=False):
    '''  Plot |phi(kx,ky)^2| as a surface plot in function of (kx,ky) and animate in time.

    The linear version is pretty useless, 2 tiles will determine the colorbar, the others remain on "0".
    '''

    if plot_start: t_range = [0,75]; time_step=0.9
    
    # Load the modules
    from stellapy.decorators.verbose_wrapper import indent  
    from stellapy.utils import ensure_dir, get_PathNameOfFolder
    import stellapy.read as read
    import numpy as np

    # Get all out.h5 files OR out.nc files inside <folder>
    netcdf_files = get_NetcdfFilesOrReducedFilesInside(folder)

    # Read the netcdf files to get potential['time', 'phi2_vs_kxky', 'phi2'] or potential[ky]['time', 'phi2_vs_kxky', 'phi2']
    potential = read.phi(folder)
    keys_shot, keys_rho = read.sorted_shot_and_rho_keys(potential, return_type='keys')

    # Read the data
    shot     = keys_shot[0]
    rho      = keys_rho[shot][0]
    ky       = list(potential[shot][rho].keys())[0]
    vec_time_full = potential[shot][rho][ky]['time']
    phi2         = potential[shot][rho][ky]['phi2_vs_kxky']
    phi2_vs_t    = potential[shot][rho][ky]['phi2']        
    dim_kx       = potential[shot][rho][ky]['dim_kx']
    vec_kx       = potential[shot][rho][ky]['vec_kx'] 
    dim_ky       = potential[shot][rho][ky]['dim_ky']
    vec_ky       = potential[shot][rho][ky]['vec_ky']  
    if (plot_phi==False) and 'time_fluxes' in potential[shot][rho][ky]:
        vec_time_flx = potential[shot][rho][ky]['time_fluxes'] 
        qflx_vs_kxky = np.abs(potential[shot][rho][ky]['qflx_kxky'])
        pflx_vs_kxky = np.abs(potential[shot][rho][ky]['pflx_kxky'])
        vflx_vs_kxky = np.abs(potential[shot][rho][ky]['vflx_kxky'])
        qflx         = np.abs(potential[shot][rho][ky]['qflx'])
        pflx         = np.abs(potential[shot][rho][ky]['pflx'])
        vflx         = np.abs(potential[shot][rho][ky]['vflx'])
    if (plot_phi==False) and not 'time_fluxes' in potential[shot][rho][ky]:
        print("EXIT: there are no fluxes found in the files.")
        return
    if plot_phi==False:
        phi2_vs_t = qflx
        phi2 = qflx_vs_kxky
        vec_time_full = vec_time_flx

    # Only plot a specific time range
    if t_range:
        phi2 = phi2[vec_time_full > t_range[0]];        vec_time = vec_time_full[vec_time_full > t_range[0]] 
        phi2 = phi2[vec_time_full < t_range[1]];        vec_time = vec_time_full[vec_time_full < t_range[1]]
        dim_time = len(vec_time)
    else: 
        vec_time = vec_time_full
        dim_time = len(vec_time)

    # Only plot one frame every time_step
    if time_step:
        vec_time_temp = [vec_time[0]]
        phi2_temp = [phi2[0, :, :]]
        t_previous = vec_time[0]
        for i in range(len(vec_time)):
            if ((vec_time[i] - t_previous) >= time_step):
                vec_time_temp.append(vec_time[i])
                phi2_temp.append(phi2[i, :, :])
                t_previous = vec_time[i]            
        vec_time = np.array(vec_time_temp)
        phi2 = np.array(phi2_temp)
        dim_time = len(vec_time)

    # Interpolate the frames
    vec_kx, vec_ky, phi2 = interpolate_y_over_kxky(phi2, vec_time, interpolate, vec_kx, vec_ky)

    # Labels for the plots
    xlabel = '$k_x \\rho_i$'
    ylabel = '$k_y \\rho_i$'
    title  = '$|\\varphi(k_x, k_y)|^2$' if plot_phi==True else '$Q/Q_{gB}$'

    # Make sure the folders for the movie exist
    directory_movie = get_PathNameOfFolder(folder) + '/' + 'movies'
    ensure_dir(directory_movie)
    if log == True and plot_phi==True:   movie_file = directory_movie + '/phi2_vs_kxky_'+cmap+'_log.mp4'
    if log == False and plot_phi==True:  movie_file = directory_movie + '/phi2_vs_kxky_'+cmap+'_lin.mp4'
    if log == True and plot_phi==False:  movie_file = directory_movie + '/q_vs_kxky_'+cmap+'_log.mp4'
    if log == False and plot_phi==False: movie_file = directory_movie + '/q_vs_kxky_'+cmap+'_lin.mp4'
    if time_step:    movie_file = movie_file.split('.mp4')[0] + '_dt' + str(time_step) + '.mp4'
    if t_range:      movie_file = movie_file.split('.mp4')[0] + '_t'+ str(t_range[0]) + '-' + str(t_range[1]) + '.mp4'

    # Make the movie
    movie_2d(phi2, phi2_vs_t, vec_kx, vec_ky, vec_time, vec_time_full, dim_time-1,plot_phi,plot_start,movie_file,xlabel,ylabel,title,cmap=cmap,log=log) 


# stellapy.plotbox.movie_2d
def movie_2d(phi2,phi2_t,xin,yin,vec_time,vec_time_full,nframes,plot_phi,plot_start,outfile,xlab='',ylab='',title='',cmap='RdBu',log=True):
    ''' Plot a surface plot of phi2 in function of (x,y) with xin and yin the vectors along the axis. '''

    # Load the modules
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import animation
    from matplotlib.colors import LogNorm
    import matplotlib.gridspec as gridspec
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from stellapy.decorators.verbose_wrapper import indent 

    # Remove broken tiles
    print(indent, 'Warning: all tiles with phi2<E-33 or phi2>100 are removed and shown in black.')
    phi2_nonan = phi2
    phi2[phi2[:,:,:]<1.E-30] = np.nan
    phi2[phi2[:,:,:]>100] = np.nan

    # Initiate the max and min of phi(kx,ky) for each time value
    phi2max = np.arange(len(vec_time),dtype=float)
    phi2min = np.arange(len(vec_time),dtype=float)

    # Find the maximum, minimum and mean of phi(kx,ky) for each time value
    for i in range(len(vec_time)):
        phi2max[i] = np.absolute(np.nanmax(phi2[i,:,:]))
        phi2min[i] = np.absolute(np.nanmin(phi2[i,:,:]))

    # Show nan values in black
    cmap = plt.get_cmap(cmap)
    cmap.set_bad(color='black')

    # Create the figure, each frame uses this figure
    fig = plt.figure(figsize=(18,6))
    fig.tight_layout(pad=0.15, h_pad=0, w_pad=0)
    gs1 = gridspec.GridSpec(1, 2)
    gs1.update(top=0.9, bottom = 0.15, wspace=0.45, left=0.07, right=0.9)
    ax1 = plt.subplot(gs1[0])
    ax2 = plt.subplot(gs1[1])
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.1)
    extent = [xin.min(),xin.max(),yin.min(),yin.max()]
    vmin = max(phi2min.mean(), 1.E-22)  
    vmax = phi2max.max()
    if plot_start: vmin = max(phi2min.min(), 1.E-22)  
    if not plot_start and not plot_phi: 
        vmin = 1.E-3; vmax = 1.E-1
    if not plot_start and plot_phi: 
        vmin = 1.E-2; vmax = 1.E1

    # Plot the first frame: the phi(kx,ky)
    if log == True:   norm = LogNorm(vmin=vmin, vmax=vmax)
    if log == False:  norm = None
    im1 = ax1.imshow(np.transpose(phi2[0,:,:]),cmap=cmap,vmin=vmin,vmax=vmax,extent=extent,\
                    interpolation='nearest',origin='lower',aspect='auto',norm=norm)  
    plt.colorbar(im1,cax=cax) 
    ax1.set_xlabel(xlab)
    ax1.set_ylabel(ylab)
    # Plot the first frame: the title
    title_text = '\ \ \ \ \ \ \ \ '+title+'\ \ \ \ t = '+"{0:3.0f}".format(vec_time[0])
    im_title = ax1.text(0.5,1.02,title_text, size=plt.rcParams["axes.titlesize"], ha="center", transform=ax1.transAxes)
    # Plot the first frame: 
    color = 'black'
    im2, = ax2.plot(vec_time_full,phi2_t,color=color)
    ax2.set_xlabel("$t\, v_{\mathrm{th}}/a$")
    if plot_phi==True: ax2.set_ylabel("$|\\varphi(t)|^2$", color=color) 
    if plot_phi==False: ax2.set_ylabel("$Q/Q_{gB}$", color=color)   
    ax2.tick_params(axis='y', labelcolor=color)
    ax2.autoscale()
    ax2.set_yscale("log")
    ax2.set_xlim(xmin=0, xmax=400)
    ax2.grid()
    im3 = ax2.axvline(x=vec_time[0], color='black')


    # Initiate an array to hold the frames
    ims = []
    ims.append([im1, im2, im3, im_title])

    # Add the other frames
    for i in range(1,nframes):
        im1 = ax1.imshow(np.transpose(phi2[i,:,:]),cmap=cmap,vmin=vmin,vmax=vmax,extent=extent,\
             interpolation='nearest',origin='lower',aspect='auto',norm=norm)
        im3 = ax2.axvline(x=vec_time[i], color='black')
        title_text = '\ \ \ \ \ \ \ \ '+title+'\ \ \ \ t = '+"{0:3.0f}".format(vec_time[i])
        im_title   = ax1.text(0.5,1.02, title_text, size=plt.rcParams["axes.titlesize"], ha="center", transform=ax1.transAxes)
        ims.append([im1, im2, im3, im_title])

    # Save the GIF
    ani = animation.ArtistAnimation(fig,ims,interval=50,blit=True)
    print("\n The animation will be saved to", outfile)
    print(" Compiling might take a long time...")
    ani.save(outfile)
    print(" The animation is saved to", outfile)



#================================================================
# Interpolate
#================================================================

def interpolate_y_over_kxky(y_vs_kxky, vec_time, interpolate, vec_kx, vec_ky):
    from stellapy.decorators.verbose_wrapper import indent  

    # Remove broken tiles
    print(indent, 'Warning: all tiles with phi2<E-30 or phi2>100 are removed and shown in black in phi'+str(np.shape(y_vs_kxky))+'.')
    print(indent, '   Tiles lower than E-30 are at:', np.argwhere(y_vs_kxky[0,:,:]<1.E-30))
    print(indent, '   Their values are:', y_vs_kxky[0, y_vs_kxky[0,:,:]<1.E-30])
    print(indent, '   Tiles bigger than 100 are at:  ', np.argwhere(y_vs_kxky[0,:,:]>100))
    print(indent, '   Their values are:', y_vs_kxky[0, y_vs_kxky[0,:,:]>100])
    filter_low  = y_vs_kxky[:,:,:]<1.E-30 # kx,ky = 0,0 is 2.58485853e-32
    filter_high = y_vs_kxky[:,:,:]>100 # kx,ky = left right,0 is 28
    
 
    # Interpolate the data
    if interpolate: 
        y_vs_kxky[filter_low] = np.nan
        y_vs_kxky[filter_high] = np.nan
        xnew = np.linspace(vec_ky[0], vec_ky[-1], int(len(vec_ky))*interpolate)
        ynew = np.linspace(vec_kx[0], vec_kx[-1], int(len(vec_kx))*interpolate)
        vec_ky_interp, vec_kx_interp = np.meshgrid(xnew, ynew)
        vec_ky_interp = vec_ky_interp[0,:]; vec_kx_interp = vec_kx_interp[:,0]
        y_vs_kxky_interp = np.zeros((len(y_vs_kxky),len(vec_kx_interp), len(vec_ky_interp)))

        for i in range(len(y_vs_kxky[:,0,0])):
            for index in np.argwhere(np.isnan(y_vs_kxky[i,:,:])): # At t,ky,kx=i,0,0 the values are nan's
                if ( index[0]<=27 or index[0]==len(y_vs_kxky[i,:,1])-1 ) and (not index[0]==0): 
                    y_vs_kxky[i, index[0], index[1]] = y_vs_kxky[i, index[0]-1, index[1]]  
                if ( index[0]>=28 or index[0]==0 ) and (not index[0]==len(y_vs_kxky[i,:,1])-1):           
                    y_vs_kxky[i, index[0], index[1]] = y_vs_kxky[i, index[0]+1, index[1]]
            function = interp2d(vec_ky, vec_kx, y_vs_kxky[i,:,:], kind='linear')
            y_vs_kxky_interp[i,:,:] = function(xnew,ynew)

    # Remove broken tiles
    if interpolate: 
        for i in range(len(y_vs_kxky[:,0,0])):
            y_vs_kxky_interp[np.repeat(np.repeat(filter_low,  interpolate, axis=2), interpolate, axis=1)] = np.nan
            y_vs_kxky_interp[np.repeat(np.repeat(filter_high, interpolate, axis=2), interpolate, axis=1)] = np.nan
    if not interpolate:
        y_vs_kxky[filter_low] = np.nan
        y_vs_kxky[filter_high] = np.nan

    if not interpolate:     return vec_kx, vec_ky, y_vs_kxky
    if interpolate:         return vec_kx_interp, vec_ky_interp, y_vs_kxky_interp

