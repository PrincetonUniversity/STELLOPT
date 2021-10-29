
    
# Load modules
import sys, os
import numpy as np
import matplotlib as mpl
from scipy.signal import find_peaks, peak_widths
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# from .plotbox import plotbox_2d
# from stellapy.utils import get_filesinfolder, initiate_nesteddict, get_fullpathoffolder, get_PathNameRuns
# from stellapy.utils import remove_cases_without_outncfile, sort_inputfiles_by_ky, get_unstable_modes, ensure_dir
from stellapy.utils.decorators import verbose_wrapper
from mpl_toolkits.axes_grid1 import make_axes_locatable

#================================================================
# Write the length of phi
#================================================================

@verbose_wrapper 
# stellapy.plot.phi
def plot_potentialVsModesVsZ(folder, kmin=0, kmax=1000, save=False):
    ''' 
    Notes
    -----
    Structure of the *.final_fields file
    # z   z-thet0   aky   akx   real(phi)   imag(phi)   real(apar)   imag(apar)   z_eqarc-thet0

    '''

#================================================================
# Get the electrostatic potential squared
#================================================================

    # Store the data for the plots in a nested dictionary data[kx][ky][quantity]
    data       = initiate_nesteddict()
    stableness = initiate_nesteddict()

    # Read all the case folders
    folders=["th_w7xr348+252_0004_nf5", "th_w7xr348+252_0004_nf10","th_w7xr348+252_0004_nf15",\
                  "th_w7xr348+252_0004_nf20","th_w7xr348+252_0004_nf25","th_w7xr348+252_0004_nf30"]
#    folders=["th_w7xr348+252_0004_lowres", "th_w7xr348+252_0004_nf15",\
#                  "th_w7xr348+252_0004_nf20","th_w7xr348+252_0004_nf25","th_w7xr348+252_0004_nf30"]
    folders=[folder]
#    folders=["th_w7xr348+252_0004_nf5", "th_w7xr348+252_0004_nf15","th_w7xr348+252_0004_nf30"]
#    folders=["th_w7xr348+252_0004_nf20", "th_w7xr348+252_0004_nf25","th_w7xr348+252_0004_nf30"]

    # Iterate over the folders
    for folder in folders:   

        # Get the input files inside  
        input_files = get_filesinfolder(folder, end='.in')

        # Remove the input files who don't have a corresponding *out.nc file
        input_files = remove_cases_without_outncfile(folder, input_files)

        # Sort input_files by ky and only keep the files fullfilling: kmin < ky < mkax
        input_files, vec_ky = sort_inputfiles_by_ky(folder, input_files, kmin, kmax)

        # Get the data for the plots
        for input_file in input_files:

            # Read all the relevant variables from the netcdf file of the case
            netcdf_data     = read.netcdf(folder, input_file)
            dim_time        = netcdf_data['dim_time']                                          
            dim_kx          = netcdf_data['dim_kx']
            dim_ky          = netcdf_data['dim_ky']
            dim_z           = netcdf_data['dim_z']
            vec_z           = netcdf_data['vec_z']
            vec_ky          = netcdf_data['vec_ky']
            vec_kx          = netcdf_data['vec_kx']

            # Read the nfield_periods
            nfield_periods  = read.input_grid(folder, input_file)['nfield_periods']

            # For comparing the same ky
            if 'test' in input_file and len(vec_ky) == 1:
                vec_ky      = list(vec_ky)
                vec_ky[0]   = str(input_file.split('_test_')[1].split('.')[0].replace('_',' '))

            # Read the *.final_fields file: [z   z-thet0  ky  kx  real(phi)  imag(phi)  real(apar) imag(apar) z_eqarc-thet0]
            fields_file  = get_filesinfolder(folder, input_file=input_file, end='final_fields')
            if len(get_filesinfolder(folder,end='.geometry')) == 1:
                geometry_file = get_filesinfolder(folder,end='.geometry')[0]
            else:
                geometry_file = get_filesinfolder(folder,input_file=input_file,end='.geometry')
            if fields_file is not None and geometry_file is not None:
                fields_file  = get_fullpathoffolder(folder) + '/' + fields_file
                try:    fields_data  = np.loadtxt(fields_file, dtype='float').reshape(dim_z, dim_kx, dim_ky, 10) # Old code has 10 columns
                except: 
                    try:    fields_data  = np.loadtxt(fields_file, dtype='float').reshape(dim_z, dim_kx, dim_ky, 11) # New code has 11 columns
                    except: 
                        try:    fields_data  = np.loadtxt(fields_file, dtype='float').reshape(513, dim_kx, dim_ky, 11) 
                        except: fields_data  = np.loadtxt(fields_file, dtype='float').reshape(513, dim_kx, dim_ky, 10) 

                # Get the stability: either plot only the unstable modes if stable = False; or only the stable modes if stable = True
                vec_ky_unstable = get_unstable_modes(folder, input_file, vec_ky, False, stableness, 0.0)

                # Read zeta
                geometry_path = get_fullpathoffolder(folder) + '/' +  geometry_file
                geometry_data = np.loadtxt(geometry_path, dtype='float', skiprows=4).reshape(-1, 13)
                for ky in vec_ky:
                    data[nfield_periods][ky]['z'] = geometry_data[:,2]

                # Store the data for each input_files in <data> and calculate the length of phi
                for kx in vec_kx:
                    for ky in vec_ky:
                        kx_i = list(vec_kx).index(kx)
                        ky_i = list(vec_ky).index(ky)
                        data[nfield_periods][ky]['ky'] = ky
                        # Read Real(phi) and Imag(phi)    
                        phi_real    = fields_data[:,kx_i,ky_i,4]
                        phi_imag    = fields_data[:,kx_i,ky_i,5]
                        phi         = [ complex(phi_real[i],phi_imag[i]) for i in range(len(phi_real)) ]
                        if   np.max(abs(np.array(phi))) > 1.E300: phi = [ phi*1.E-300 for phi in phi] # Otherwise the square is infinite
                        elif np.max(abs(np.array(phi))) > 1.E200: phi = [ phi*1.E-200 for phi in phi] # Otherwise the square is infinite
                        elif np.max(abs(np.array(phi))) > 1.E100: phi = [ phi*1.E-100 for phi in phi] # Otherwise the square is infinite
                        phi2        = [ abs(phi*phi.conjugate()) for phi in phi ]
                        # Store it
                        data[nfield_periods][ky]['phi2']        = phi2
                        data[nfield_periods][ky]['phi']         = phi_real

#================================================================
# Get the magnetic field and curvature
#================================================================

    # Get the geometry files inside  
    geometry_files = get_filesinfolder(folder, end='.geometry')[::-1]    
    geometry_files = [ geometry_file for geometry_file in geometry_files if int(geometry_file.split('ky_')[1].split('.')[0]) > 10 ] 
    # Store the data for the plots in a nested dictionary
    keys  = ["alpha","zed", "zeta", "bmag", "gradpar", "gds2", "gds21", "gds22", "gds23", "gds24", "gbdrift", "cvdrift", "gbdrift0"]
    data2 = initiate_nesteddict()

    # Get the data for the plots
    for geometry_file in geometry_files:

        # Read all the relevant variables from the geometry file:
        # alpha; zed; zeta; bmag; gradpar; gds2; gds21; gds22; gds23; gds24; gbdrift; cvdrift; gbdrift0
        geometry_path = get_fullpathoffolder(folder) + '/' +  geometry_file
        geometry_data = np.loadtxt(geometry_path, dtype='float', skiprows=4).reshape(-1, 13)

        # Read the magnetic field identifier
        magnetic_field = str(geometry_file.split('B')[1].split('_')[0])

        # Store the data
        for key in keys: 
            data2[magnetic_field][key] = geometry_data[:,keys.index(key)]
        

#================================================================
# plot_length_phi
#================================================================
   
    # Labels for Re[phi](t) and Im[phi](t)
    x_label = '$\\zeta$'
    y_label = '$k_{y}\\rho_i$' 

    # Create the figure
    if len(folders)==6:
        fig = plt.figure(figsize=(20,25))
        gs1 = gridspec.GridSpec(6, 1)
        gs1.update(wspace=0.025, hspace=0.00, top=0.92, right=0.7, bottom = 0.1)
        ax1 = fig.add_subplot(gs1[0])
        ax2 = fig.add_subplot(gs1[1])
        ax3 = fig.add_subplot(gs1[2])
        ax4 = fig.add_subplot(gs1[3])
        ax5 = fig.add_subplot(gs1[4])
        ax6 = fig.add_subplot(gs1[5])
        plotbox_2d(x_label=x_label, y_label=y_label, ax=ax1) 
        plotbox_2d(x_label=x_label, y_label=y_label, ax=ax2) 
        plotbox_2d(x_label=x_label, y_label=y_label, ax=ax3) 
        plotbox_2d(x_label=x_label, y_label=y_label, ax=ax4) 
        plotbox_2d(x_label=x_label, y_label=y_label, ax=ax5) 
        plotbox_2d(x_label=x_label, y_label=y_label, ax=ax6) 

    if len(folders)==3:
        fig = plt.figure(figsize=(20,7))
        gs1 = gridspec.GridSpec(2, 100)
        gs1.update(wspace=0.025, hspace=0.2, top=0.92, right=0.7, bottom = 0.15)
        ax1 = fig.add_subplot(gs1[0,:25])
        ax2 = fig.add_subplot(gs1[0,26:])
        ax3 = fig.add_subplot(gs1[1,:])
        plotbox_2d(x_label=x_label, y_label=y_label, ax=ax1) 
        plotbox_2d(x_label=x_label, y_label="", ax=ax2) 
        plotbox_2d(x_label=x_label, y_label=y_label, ax=ax3) 

    if len(folders)==1:
        plt.rc('font', size=20)
        figsize=(6,4)
        figsize=(10,4)
        fig = plt.figure(figsize=figsize)
        gs1 = gridspec.GridSpec(100, 1)
        gs1.update(wspace=0.025, hspace=0.05, top=0.92, right=0.7, bottom = 0.15)
        ax1 = fig.add_subplot(gs1[0:79])
        plotbox_2d(x_label="", y_label=y_label, ax=ax1) 
        ax2 = fig.add_subplot(gs1[80:89])
        plotbox_2d(x_label="", y_label="", ax=ax2) 
        ax3 = fig.add_subplot(gs1[90:100])
        plotbox_2d(x_label="", y_label="", ax=ax3) 
    

    # One color for each nfield_period
    nfield_periods_values = list(data.keys())
    if len(folders)==6: 
        nfield_periods_values = np.array(list(data.keys()))
        nfield_periods_values = list(nfield_periods_values[np.argsort(nfield_periods_values)[::-1]])
    color_values = plt.cm.jet( np.linspace( 0,1,len(nfield_periods_values) ) ) 
    color_values = plt.cm.inferno( np.linspace( 0,1,101 ) ) [::-1]
    marker_styles = ['o','X','D','P','s','x','d','p']*5 

    # Plot Re[phi](t) and Im[phi](t) for the chosen kx value and either the stable or not-stable modes
    if len(folders) > 1:
        for nfield_periods in nfield_periods_values:
            if nfield_periods_values.index(nfield_periods) == 0: ax=ax1; print("ax1: ", nfield_periods)
            if nfield_periods_values.index(nfield_periods) == 1: ax=ax2; print("ax2: ", nfield_periods)
            if nfield_periods_values.index(nfield_periods) == 2: ax=ax3; print("ax3: ", nfield_periods)
            if nfield_periods_values.index(nfield_periods) == 3: ax=ax4; print("ax4: ", nfield_periods)
            if nfield_periods_values.index(nfield_periods) == 4: ax=ax5; print("ax5: ", nfield_periods)
            if nfield_periods_values.index(nfield_periods) == 5: ax=ax6; print("ax6: ", nfield_periods)

        # Create the coordinate map
        ky_values = list(data[nfield_periods].keys())
        first_ky = list(data[nfield_periods].keys())[0]
        z_radials = data[nfield_periods][first_ky]['z'] 
        height = np.zeros((len(z_radials), len(ky_values)))

        for ky in data[nfield_periods].keys():
            # Unpack the data
            phi2 = data[nfield_periods][ky]['phi2']  
            z_radials = data[nfield_periods][ky]['z'] 
            height[:,ky_values.index(ky)] = phi2/np.max(phi2)

        ax.pcolormesh(z_radials, ky_values, np.transpose(height), cmap=plt.get_cmap('inferno_r'), norm=None)   

    if len(folders) == 1:
        ax = ax1

        for nfield_periods in nfield_periods_values:
            # Create the coordinate map
            ky_values = list(data[nfield_periods].keys())
            if nfield_periods==11.16976 and 4.5 in ky_values: ky_values.append(5.0)
            first_ky = list(data[nfield_periods].keys())[0]
            z_radials = data[nfield_periods][first_ky]['z'] 
            height = np.zeros((len(z_radials), len(ky_values)))
            print(nfield_periods, ky_values)

            for ky in data[nfield_periods].keys():
                # Unpack the data
                phi2 = data[nfield_periods][ky]['phi2']  
                z_radials = data[nfield_periods][ky]['z'] 
                height[:,ky_values.index(ky)] = phi2/np.max(phi2) 

            ax.pcolormesh(z_radials, ky_values, np.transpose(height), cmap=plt.get_cmap('inferno_r'), norm=None)   

        # Set up the axis
        ax.set_ylim(ymin=0.0, ymax=20)
        #ax.set_ylim(ymin=0.0, ymax=5)
        ax.xaxis.set_major_locator(plt.MultipleLocator(np.pi / 1))
        ax.xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 4))
        ax.xaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
        ax.set_xlim(xmin=-6.6*np.pi, xmax=6.6*np.pi)
        if len(folders)==3:
            if nfield_periods > 5 and len(folders)==3:
                ax.xaxis.set_major_locator(plt.MultipleLocator(np.pi / 1))
                ax.xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 4))
                ax.xaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
                ax.set_xlim(xmin=-1.1*np.pi, xmax=1.1*np.pi)
            if nfield_periods > 10 and len(folders)==3:
                ax.xaxis.set_major_locator(plt.MultipleLocator(np.pi / 1))
                ax.xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 4))
                ax.xaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
                ax.set_xlim(xmin=-2.2*np.pi, xmax=2.2*np.pi)
            if nfield_periods > 15 and (len(folders)==3 or len(folders)==1):
                ax.xaxis.set_major_locator(plt.MultipleLocator(np.pi / 1))
                ax.xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 4))
                ax.xaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
                ax.set_xlim(xmin=-3.4*np.pi, xmax=3.4*np.pi)
            if nfield_periods > 30 and len(folders)==3:
                ax.xaxis.set_major_locator(plt.MultipleLocator(np.pi / 1))
                ax.xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 4))
                ax.xaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
                ax.set_xlim(xmin=-6.6*np.pi, xmax=6.6*np.pi)
        if len(folders)==1:
            ax.xaxis.set_major_locator(plt.MultipleLocator(np.pi / 1))
            ax.xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 4))
            ax.xaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
            ax.set_xlim(xmin=-3.4*np.pi, xmax=3.4*np.pi)

        # Finish the figures
        if len(folders)!=3 or (len(folders)==3 and nfield_periods>12):
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size='5%', pad=0.1) 
            cmap = mpl.cm.get_cmap('inferno_r')
            norm = mpl.colors.Normalize(vmin=0, vmax=1)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=1)) 
            sm._A = []
            cbar = plt.colorbar(sm, cax=cax)
            cbar.set_ticks([])
            #cbar.ax.set_yticklabels(["{:.0%}".format(i) for i in cbar.get_ticks()]) # set ticks of your format %.0f%%'
            cbar.set_label('$\\varphi^2/\\varphi^2_{max}$', rotation=270, labelpad=40)
        if len(folders)==3 and (nfield_periods>12 and nfield_periods<30):
            ax.yaxis.set_major_formatter(plt.NullFormatter())
    #ax1.legend(loc='lower right', shadow=True, ncol=4, labelspacing=0.0, prop={'size':24}, )
    #ax1.legend(loc='upper center', bbox_to_anchor=(1, -0.15), labelspacing=0.0, ncol=5, shadow=True, prop={'size':24})

    # Add the shape of the magnetic field and the curvature
    if len(folders)==1:
        for magnetic_field in data2.keys():
            for key in keys:
                maximum1 = np.max(np.abs(data2[magnetic_field]['cvdrift']))
                minimum3 = np.min(np.abs(data2[magnetic_field]['bmag']))
                maximum3 = np.max(np.abs(data2[magnetic_field]['bmag']-minimum3))
                if   magnetic_field == "348": sign=1
                elif magnetic_field == "169": sign=-1
                else:                         sign=1

                # Create the coordinate map
                ky_values = [0,1]
                z_radials = data2[magnetic_field]['zeta']
                gbdrift = np.zeros((len(z_radials), len(ky_values)))
                bmag    = np.zeros((len(z_radials), len(ky_values)))
                for ky in ky_values: 
                    gbdrift[:,ky_values.index(ky)]  = data2[magnetic_field]['gbdrift']/maximum1*sign
                    bmag[:,ky_values.index(ky)]     = (data2[magnetic_field]['bmag']-minimum3)/maximum3 
                ax2.pcolormesh(z_radials, ky_values, np.transpose(bmag),    cmap=plt.get_cmap('inferno_r'), vmin=0, vmax=1)  
                ax3.pcolormesh(z_radials, ky_values, np.transpose(gbdrift), cmap=plt.get_cmap('bwr_r'), norm=None, vmin=-1, vmax=1)  
    
        for ax in [ax2,ax3]:
            ax.set_ylim(ymin=0, ymax=1)
            ax.xaxis.set_major_locator(plt.MultipleLocator(np.pi / 1))
            ax.xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 4))
            ax.xaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
            ax.set_xlim(xmin=-3.4*np.pi, xmax=3.4*np.pi)
            ax.yaxis.set_minor_locator(ticker.FixedLocator(np.linspace(-2,2)))
            ax.yaxis.set_major_locator(ticker.FixedLocator([-2,2]))
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size='5%', pad=0.1) 
            if ax == ax2:   cmap = mpl.cm.get_cmap('inferno_r')
            if ax == ax3:   cmap = mpl.cm.get_cmap('bwr_r')
            if ax == ax2:   norm = plt.Normalize(vmin=0, vmax=1)
            if ax == ax3:   norm = plt.Normalize(vmin=-1, vmax=1)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm) 
            sm._A = []
            cbar = plt.colorbar(sm, cax=cax)
            cbar.set_ticks([])
            if ax == ax3:   cbar.set_label("$B \ \ \\mathcal{K}_2$", rotation=270, labelpad=40)
            if ax == ax2:   cbar.set_label("", rotation=270, labelpad=40)
        
        for ax in [ax1, ax2, ax3]:
            ax.set_xlim(xmin=-3.45*np.pi, xmax=3.45*np.pi)
            ax.set_xlim(xmin=-2.25*np.pi, xmax=2.25*np.pi)
        for ax in [ax1, ax2]:
            ax.xaxis.set_major_formatter(plt.NullFormatter())

    #============================
    # Save the plots to a folder 
    #============================
    # Load the modules
    from stellapy.decorators.verbose_wrapper import indent

    if save and (os.getcwd().split('/')[1] != 'marconi_work'):
        # Make sure the folders for the figures exist
        directory_png = get_fullpathoffolder(folder) + '/' + 'figures_png'
        directory_eps = get_fullpathoffolder(folder) + '/' + 'figures_eps'
        ensure_dir(directory_eps);  ensure_dir(directory_png);
        # Get the name of the figure
        figure_name = "Shape_phi_below5"
        print("\n", indent, "The plot is saved as ", figure_name, sep="")
        #save the figure
        #plt.savefig(directory_eps+'/'+figure_name+'.eps', format='eps')
        plt.savefig(directory_png+'/'+figure_name+'.png', format='png')
    plt.show()
    return


#================================================================
# formatter
#================================================================

def multiple_formatter(denominator=2, number=np.pi, latex='\pi'):
    def gcd(a, b):
        while b:
            a, b = b, a%b
        return a
    def _multiple_formatter(x, pos):
        den = denominator
        num = np.int(np.rint(den*x/number))
        com = gcd(num,den)
        (num,den) = (int(num/com),int(den/com))
        if den==1:
            if num==0:
                return r'$0$'
            if num==1:
                return r'$%s$'%latex
            elif num==-1:
                return r'$-%s$'%latex
            else:
                return r'$%s%s$'%(num,latex)
        else:
            if num==1:
                return r'${%s}/{%s}$'%(latex,den)
            elif num==-1:
                return r'${-%s}/{%s}$'%(latex,den)
            else:
                return r'${%s%s}/{%s}$'%(num,latex,den)
    return _multiple_formatter


class Multiple:
    def __init__(self, denominator=2, number=np.pi, latex='\pi'):
        self.denominator = denominator
        self.number = number
        self.latex = latex

    def locator(self):
        return plt.MultipleLocator(self.number / self.denominator)

    def formatter(self):
        return plt.FuncFormatter(multiple_formatter(self.denominator, self.number, self.latex))



#            peaks = data[nfield_periods][ky]['peaks']       
#            phi_half = data[nfield_periods][ky]['phi_half']    
#            start_half = data[nfield_periods][ky]['start_half']   
#            end_half = data[nfield_periods][ky]['end_half']     
#            phi_full = data[nfield_periods][ky]['phi_full']     
#            start_full = data[nfield_periods][ky]['start_full']  
#            end_full = data[nfield_periods][ky]['end_full']      
            # Plot it
#            for i in range(len(peaks)): # alpha=(phi_real/np.max(phi_real))[peaks[i]]
#                #ax.plot(z_radials[peaks[i]], ky, "x", color=color_values[int((phi_real/np.max(phi_real))[peaks[i]]*100)])
#                ax.hlines(ky, start_half[i], end_half[i], color=color_values[int((phi_real/np.max(phi_real))[peaks[i]]*100)], lw=lw)
#                #ax.hlines(ky, start_full[i], end_full[i], color=color_values[int((phi_real/np.max(phi_real))[peaks[i]]*100)], lw=5)
#                #ax.hlines([ky]*len(start_full), start_full, end_full, color=color)
#            if False: 
#                phi_norm = int(phi_real/np.max(phi_real)*100)
#                color_line = [color_values[phi_norm[i]] for i in len(phi_norm)]
#                scatter([ky]*len(phi_real),phi_real,color_line,'fill')

