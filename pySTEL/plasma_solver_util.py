#!/usr/bin/env python3
# -*- coding: utf-8 -*-


if __name__=="__main__":
	import sys
	from argparse import ArgumentParser
	import matplotlib.pyplot as plt
	import matplotlib
	import numpy as np
	import joblib
	from datetime import datetime
	from libstell.fusion import FUSION,EC
	#from libstell.POPCON import POPCON
	from scipy.integrate import cumulative_trapezoid
	parser = ArgumentParser(description= 
		'''Utility for plotting 1D transport simulations.''')
	parser.add_argument("--output", dest="output_ext",
		help="Output joblib file name.", default = None)
	parser.add_argument("-p", "--plot", dest="lplot", action='store_true',
		help="Make plots.", default = False)
	parser.add_argument("--plot_popcon", dest="lplot_popcon", action='store_true',
		help="Overplot run on POPCON plot.", default = False)
	parser.add_argument("--save", dest="lsave", action='store_true',
		help="Save the plots with ext names.", default = False)
	args = parser.parse_args()
	if args.output_ext: 
		solver = joblib.load(f'{args.output_ext}.joblib')
		try:     
			Area = solver.dVdr(solver.rho_grid)
		except:
			Area = solver.dVdr
		try:
		    S_ECRH = solver.explicit_energy_sources['electrons']['external_gaussian'][:,:]
		except:
		    S_ECRH = solver.explicit_energy_sources['electrons']['time_dependent_gaussian'][:,:]
		S_alpha = 0.0
		for species in solver.list_of_species:
		    try:
		        S_alpha += solver.explicit_energy_sources[species]['alpha_heating'][:,:]
		    except:
		        S_alpha += 0       
		S_Bremm = solver.explicit_energy_sources['electrons']['Bremsstrahlung'][:,:]
		P_ECRH  = cumulative_trapezoid(S_ECRH*Area, solver.r_grid,axis=1,initial=0.0)
		P_alpha = cumulative_trapezoid(S_alpha*Area,solver.r_grid,axis=1,initial=0.0)
		P_Bremm = cumulative_trapezoid(S_Bremm*Area,solver.r_grid,axis=1,initial=0.0)
		P_TOTAL = P_ECRH+P_alpha
		try:
		    S_fueling = solver.explicit_particle_sources['deuterium']['external_gaussian'][:,:]
		except:
		    #S_fueling = solver.explicit_particle_sources['deuterium']['time_dependent_gaussian'][:,:]
		    S_fueling = solver.explicit_particle_sources['deuterium']['PID_edense_gaussian'][:,:]
		pellet_content_times_freq = 5E20/0.1
		P_fueling = np.trapezoid(S_fueling*Area, solver.r_grid, axis=1)
		fusion  = FUSION()
		tauiss04 = lambda a,R,P,n_avg,B,iota: 0.134 * a**2.28 * R**0.64 * (P/1e6)**-0.61 * (n_avg/1e19)**0.54 * B**0.84 * iota**0.41
		aminor = solver.aminor
		Rmajor = solver.Rmajor
		try:
		    B = solver.B[:,0] # use axes value
		except:
		    B = solver.B
		try:
		    ir2o3 = np.argmin(np.abs(solver.rho_grid-2/3))
		    iota = solver.iota[:,ir2o3]
		except:
		    print('WARNING: Using iota=0.9 to compute ISS04')
		    iota = 0.9
		n_avg = {}
		for species in solver.list_of_species:
		    n_avg[species] = np.trapezoid(solver.N[species][:,:]*Area,solver.r_grid,axis=1) / np.trapz(Area,solver.r_grid)
		try:
		    tau_ISS04 = [tauiss04(aminor[it],Rmajor[it],P_ECRH[it,-1]+P_alpha[it,-1],n_avg['electrons'][it],B[it],iota[it]) for it,_ in enumerate(solver.time)]
		except:
		    tau_ISS04 = [tauiss04(aminor,Rmajor,P_ECRH[it,-1]+P_alpha[it,-1],n_avg['electrons'][it],B,iota) for it,_ in enumerate(solver.time)]
		pressure = 0.0
		for species in solver.list_of_species:
		    pressure += 1.5*solver.N[species][:,:]*solver.T[species][:,:]*EC
		    #pressure += (3.0/2.0)*solver.N[species][:,:]*solver.T[species][:,:]*EC
		W_total = np.trapezoid(pressure*Area,solver.r_grid,axis=1)
		dWdt = np.gradient(W_total,solver.time)
		tau_E = W_total / (-dWdt+P_ECRH[:,-1]+P_alpha[:,-1])
		iss04_fact = tau_E[-1]/tau_ISS04[-1]
		if args.lplot:
			px = 1/plt.rcParams['figure.dpi']
			font = {'family' : 'Arial',
			        'weight' : 'normal',
			        'size'   : 18}
			matplotlib.rc('font', **font)
			#fig,ax = plt.subplots(4,1,figsize=(1800*px,2400*px))
			fig,ax = plt.subplots(4,1,figsize=(900*px,1200*px))
			ax[0].plot(solver.time,P_TOTAL[:,-1]/1E6,linewidth=2.0,color='#5faf30',label=r'$P_{\mathrm{TOTAL}}$')
			ax[0].plot(solver.time,10*P_ECRH[:,-1]/1E6,linewidth=2.0,color='blue',label=r'$P_{\mathrm{ECRH}}x10$')
			ax[0].plot(solver.time,P_alpha[:,-1]/1E6,linewidth=2.0,color='green',label=r'$P_{\mathrm{\alpha}}$')
			ax[0].plot(solver.time,-P_Bremm[:,-1]/1E6,linewidth=2.0,color='red',label=r'$P_{\mathrm{Brem.}}$')
			ax[0].grid()
			ax[0].legend()
			ax[0].set_title('GIGA')
			ax[0].set_ylabel('P [MW]')
			ax[1].plot(solver.time,solver.N['electrons'][:,0]/1E19,linewidth=2.0,color='#5faf30',label=r'$n_e$')
			ax[1].plot(solver.time,solver.N['deuterium'][:,0]/1E19,':',linewidth=2.0,color='red',label=r'$n_D$')
			ax[1].plot(solver.time,solver.N['tritium'][:,0]/1E19,linewidth=2.0,color='#004817',label=r'$n_T$')
			if 'helium4' in solver.N.keys():
			    ax[1].plot(solver.time,solver.N['helium4'][:,0]/1E19,linewidth=2.0,color='#004817',label=r'$n_{He4}$')
			if 'alphas_fast' in solver.N.keys():
			    ax[1].plot(solver.time,solver.N['alphas_fast'][:,0]/1E19,linewidth=2.0,color='green',label=r'$n_{He4-fast}$')
			if 'neon' in solver.N.keys():
			    ax[1].plot(solver.time,solver.N['neon'][:,0]/1E19,linewidth=2.0,color='magenta',label=r'$n_{Ne}$')
			ax[1].grid()
			ax[1].legend()
			ax12=ax[1].twinx()
			ax12.plot(solver.time,P_fueling/1E22,linewidth=1.0,color='black',label=r'$N$')
			ax[1].set_ylabel(r'$n_0~[10^{19}]~m^{-3}$')
			ax12.set_ylabel(r'$\dot{N}~[10^{22}]~part/s$')
			ax[2].plot(solver.time,solver.T['electrons'][:,0]/1E3,linewidth=2.0,color='#5faf30',label=r'$T_e$')
			ax[2].plot(solver.time,solver.T['deuterium'][:,0]/1E3,':',linewidth=2.0,color='red',label=r'$T_D$')
			ax[2].plot(solver.time,solver.T['tritium'][:,0]/1E3,linewidth=2.0,color='#004817',label=r'$T_T$')
			if 'helium4' in solver.N.keys():
			    ax[2].plot(solver.time,solver.T['helium4'][:,0]/1E3,linewidth=2.0,color='#004817',label=r'$T_{He4}$')
			if 'neon' in solver.N.keys():
			    ax[2].plot(solver.time,solver.T['neon'][:,0]/1E3,linewidth=2.0,color='magenta',label=r'$T_{Ne}$')
			ax[2].set_ylabel(r'$T_0~[keV]$')
			ax[2].legend()
			ax[2].grid()
			ax[3].plot(solver.time,W_total/1E9,linewidth=2.0,color='#5faf30',label=r'$W_{therm} [GJ]$')
			ax12=ax[3].twinx()
			ax12.plot(solver.time,tau_E/tau_ISS04,label=r'$\tau_{\mathrm{ISS04}}$',linewidth=2.0)
			ax12.set_ylim(0.5,1.5)
			ax[3].set_ylabel(r'$W_{therm} [GJ]$')
			ax[3].set_xlabel('Time [s]')
			ax[3].grid()
			plt.tight_layout()
			plt.show()
			if (args.lsave): fig.savefig(f'overview_{args.output_ext}.png', dpi=fig.dpi)
		if args.lplot_popcon:
			fig,ax = plt.subplots(1,1,figsize=(900*px,1200*px))
			ax.plot(solver.T['electrons'][:,0]/1E3,solver.N['electrons'][:,0]/1E20)
			ax.grid()
			ax.set_xlabel(r'$T_{e0}~[keV]$')
			ax.set_ylabel(r'$n_{e0}~\times10^{20}~[m^{-3}]$')
			ax.set_xlim(5,30)
			ax.set_ylim(0,3.0)
			plt.tight_layout()
			plt.show()



	sys.exit(0)
