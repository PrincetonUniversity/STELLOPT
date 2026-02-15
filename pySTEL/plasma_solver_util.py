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
	from libstell.popcon import POPCON
	from libstell.plasma import PLASMA
	from libstell.plasma_solver import merge_output_files
	from scipy.integrate import cumulative_trapezoid
	parser = ArgumentParser(description= 
		'''Utility for plotting 1D transport simulations.''')
	parser.add_argument("--output", dest="output_files", nargs='+',
		help="One or more output joblib files (without .joblib extension).", default = None)
	parser.add_argument("-p", "--plot", dest="lplot", action='store_true',
		help="Make plots.", default = False)
	parser.add_argument("--plot_popcon", dest="lplot_popcon", action='store_true',
		help="Overplot run on POPCON plot.", default = False)
	parser.add_argument("--plot_profs", dest="tslice_profs",
		help="Plot profiles at given timeslice.", default = -1.0, type=float)
	parser.add_argument("--save", dest="lsave", action='store_true',
		help="Save the plots with ext names.", default = False)
	args = parser.parse_args()
	if args.output_files is not None:
		output_files = [f"{f}.joblib" for f in args.output_files]
		if len(output_files) == 1:
			solver = joblib.load(output_files[0])
		else:
			solver = merge_output_files(*output_files)
		nr     = len(solver.r_grid)
		nt     = solver.Nt
		Area = solver.dVdr
		#
		keys = solver.explicit_energy_sources['electrons'].keys()
		S_ECRH = np.zeros((nt,nr))
		for ttype in ['external_gaussian','time_dependent_gaussian','PID_etemp_gaussian','PID_itemp_gaussian','PID_pfuse_gaussian']:
			if ttype in keys:
				S_ECRH += solver.explicit_energy_sources['electrons'][ttype][:,:]
		S_alpha = np.zeros((nt,nr))
		for species in solver.list_of_species:
			keys = solver.explicit_energy_sources[species].keys()
			if 'alpha_heating' in keys:
				S_alpha += solver.explicit_energy_sources[species]['alpha_heating'][:,:]
		S_Bremm = solver.explicit_energy_sources['electrons']['Bremsstrahlung'][:,:]
		P_ECRH  = cumulative_trapezoid(S_ECRH*Area, solver.r_grid,axis=1,initial=0.0)
		P_alpha = cumulative_trapezoid(S_alpha*Area,solver.r_grid,axis=1,initial=0.0)
		P_Bremm = cumulative_trapezoid(S_Bremm*Area,solver.r_grid,axis=1,initial=0.0)
		P_TOTAL = P_ECRH+P_alpha
		keys = solver.explicit_particle_sources['deuterium'].keys()
		S_fueling = np.zeros((nt,nr)) #None
		for ttype in ['external_gaussian','time_dependent_gaussian','PID_edense_gaussian','PID_pfuse_gaussian']:
			if ttype in keys:
				S_fueling += solver.explicit_particle_sources['deuterium'][ttype][:,:]
				# if type(S_fueling) is type(None):
				# 	S_fueling = solver.explicit_particle_sources['deuterium'][ttype][:,:]
				# else:
				# 	S_fueling = S_fueling + solver.explicit_particle_sources['deuterium'][ttype][:,:]
		P_fueling = np.trapezoid(S_fueling*Area, solver.r_grid, axis=1)
		fusion  = FUSION()
		tauiss04 = lambda a,R,P,n_avg,B,iota: 0.134 * a**2.28 * R**0.64 * (P/1e6)**-0.61 * (n_avg/1e19)**0.54 * B**0.84 * iota**0.41
		aminor = solver.aminor
		Rmajor = solver.Rmajor
		try:
			B = solver.B[:,0] # use axes value
		except:
			B = solver.Baxis
		try:
			ir2o3 = np.argmin(np.abs(solver.rho_grid-2/3))
			iota = solver.iota[:,ir2o3]
		except:
			iota = solver.iota23
		n_avg = {}
		for species in solver.list_of_species:
			n_avg[species] = np.trapezoid(solver.N[species][:,:]*Area,solver.r_grid,axis=1) / np.trapezoid(Area,solver.r_grid)
		try:
			tau_ISS04 = [tauiss04(aminor[it],Rmajor[it],P_ECRH[it,-1]+P_alpha[it,-1],n_avg['electrons'][it],B[it],iota[it]) for it,_ in enumerate(solver.time)]
		except:
			tau_ISS04 = [tauiss04(aminor,Rmajor,P_ECRH[it,-1]+P_alpha[it,-1],n_avg['electrons'][it],B,iota) for it,_ in enumerate(solver.time)]
		pressure = 0.0
		for species in solver.list_of_species:
			pressure += 1.5*solver.N[species][:,:]*solver.T[species][:,:]*EC
		W_total = np.trapezoid(pressure*Area,solver.r_grid,axis=1)
		dWdt = np.gradient(W_total,solver.time)
		tau_E = W_total / (-dWdt+P_ECRH[:,-1]+P_alpha[:,-1])
		iss04_fact = tau_E[-1]/tau_ISS04[-1]
		if args.lplot:
			px = 1/plt.rcParams['figure.dpi']
			font = {'family' : 'Arial',
					'weight' : 'normal',
					'size'   : 14}
			matplotlib.rc('font', **font)
			#fig,ax = plt.subplots(4,1,figsize=(1800*px,2400*px))
			fig,ax = plt.subplots(4,1,figsize=(900*px,1200*px))
			ax[0].plot(solver.time,P_TOTAL[:,-1]/1E6,linewidth=2.0,color='#5faf30',label=r'$P_{\mathrm{TOTAL}}$')
			# ax[0].plot(solver.time,10*P_ECRH[:,-1]/1E6,linewidth=2.0,color='blue',label=r'$P_{\mathrm{ECRH}}x10$')
			ax[0].plot(solver.time,P_ECRH[:,-1]/1E6,linewidth=2.0,color='blue',label=r'$P_{\mathrm{ECRH}}$')
			ax[0].plot(solver.time,P_alpha[:,-1]/1E6,linewidth=2.0,color='green',label=r'$P_{\mathrm{\alpha}}$')
			# ax[0].plot(solver.time,5*P_alpha[:,-1]/1E6,linewidth=2.0,color='k',label=r'$P_{\mathrm{fusion}}$')
			ax[0].plot(solver.time,-P_Bremm[:,-1]/1E6,linewidth=2.0,color='red',label=r'$P_{\mathrm{Brem.}}$')
			ax[0].grid()
			ax[0].legend()
			ax[0].set_ylabel('P [MW]')
			Ne = np.trapezoid(solver.N['electrons'][:,:]*Area, solver.r_grid)
			ax[1].plot(solver.time,solver.N['electrons'][:,0]/1E19,linewidth=2.0,color='#5faf30',label=r'$n_e$')
			if 'hydrogen' in solver.N.keys():
				ax[1].plot(solver.time,solver.N['hydrogen'][:,0]/1E19,':',linewidth=2.0,color='red',label=r'$n_H$')
			if 'deuterium' in solver.N.keys():
				ax[1].plot(solver.time,solver.N['deuterium'][:,0]/1E19,':',linewidth=2.0,color='red',label=r'$n_D$')
			if 'tritium' in solver.N.keys():
				ax[1].plot(solver.time,solver.N['tritium'][:,0]/1E19,linewidth=2.0,color='#004817',label=r'$n_T$')
			if 'helium4' in solver.N.keys():
				# compute fraction of helium4
				NHe4 = np.trapezoid(solver.N['helium4'][:,:]*Area, solver.r_grid)
				frac = NHe4/Ne
				ax[1].plot(solver.time,solver.N['helium4'][:,0]/1E19,linewidth=2.0,color='#004817',label=fr'$n_{{He4}}$ (f={frac[-1]*100:.2f}%)')
			if 'alphas_fast' in solver.N.keys():
				ax[1].plot(solver.time,solver.N['alphas_fast'][:,0]/1E19,linewidth=2.0,color='green',label=r'$n_{He4-fast}$')
			if 'neon' in solver.N.keys():
				# compute fraction of neon
				NNe = np.trapezoid(solver.N['neon'][:,:]*Area, solver.r_grid)
				frac = NNe/Ne
				ax[1].plot(solver.time,solver.N['neon'][:,0]/1E19,linewidth=2.0,color='magenta',label=fr'$n_{{Ne}}$ (f={frac[-1]*100:.2f}%)')
			ax[1].grid()
			ax[1].legend()
			ax12=ax[1].twinx()
			ax12.plot(solver.time,P_fueling/1E22,linewidth=1.0,color='black',label=r'$N$')
			ax[1].set_ylabel(r'$n_0~[10^{19}~m^{-3}]$')
			ax12.set_ylabel(r'$\dot{N}~[10^{22}~part/s]$')
			ax[2].plot(solver.time,solver.T['electrons'][:,0]/1E3,linewidth=2.0,color='#5faf30',label=r'$T_e$')
			if 'hydrogen' in solver.T.keys():
				ax[2].plot(solver.time,solver.T['hydrogen'][:,0]/1E3,':',linewidth=2.0,color='red',label=r'$T_H$')
			if 'deuterium' in solver.T.keys():
				ax[2].plot(solver.time,solver.T['deuterium'][:,0]/1E3,':',linewidth=2.0,color='red',label=r'$T_D$')
			if 'tritium' in solver.T.keys():
				ax[2].plot(solver.time,solver.T['tritium'][:,0]/1E3,linewidth=2.0,color='#004817',label=r'$T_T$')
			if 'helium4' in solver.T.keys():
				ax[2].plot(solver.time,solver.T['helium4'][:,0]/1E3,linewidth=2.0,color='#004817',label=r'$T_{He4}$')
			if 'neon' in solver.T.keys():
				ax[2].plot(solver.time,solver.T['neon'][:,0]/1E3,linewidth=2.0,color='magenta',label=r'$T_{Ne}$')
			ax[2].set_ylabel(r'$T_0~[keV]$')
			ax[2].legend()
			ax[2].grid()
			ax[3].plot(solver.time,W_total/1E9,linewidth=2.0,color='#5faf30',label=r'$W_{therm}$')
			ax12=ax[3].twinx()
			ax12.plot(solver.time,tau_E/tau_ISS04,label=r'$\tau_E/\tau_{\mathrm{ISS04}}$',linewidth=2.0)
			ax12.set_ylim(0.5,1.5)
			ax[3].set_ylabel(r'$W_{therm}~[GJ]$')
			ax[3].set_xlabel('Time [s]')
			ax[3].grid()
			ax[3].legend(loc='upper left')
			ax12.legend(loc='upper right')
			ax12.set_ylabel(r'$\tau_E/\tau_{\mathrm{ISS04}}$')
			plt.tight_layout()
			plt.show()
			if (args.lsave): fig.savefig(f'overview_{args.output_ext}.png', dpi=fig.dpi)
		if args.lplot_popcon:
			te_min = 2.0E3; te_max = 30.0E3
			ne_min = 1.0E19; ne_max = 3.0E20
			nte = 32; nne=32
			plasma=[[0 for x in range(nte)] for y in range(nne)] 
			te_vec = np.linspace(te_min,te_max,nte)
			ne_vec = np.linspace(ne_min,ne_max,nne)
			for i,ttemp in enumerate(te_vec):
				for j,ntemp in enumerate(ne_vec):
					plasma[i][j] = PLASMA(solver.list_of_species)
					ne0 = solver.N['electrons'][0,0]
					neE = solver.N['electrons'][0,-1]
					for k,spec in np.ndenumerate(solver.list_of_species):
						n0=ntemp*solver.N[spec][0,0]/ne0
						ne=ntemp*solver.N[spec][0,-1]/neE
						t0=ttemp
						te=t0*0.01
						# plasma[i][j].set_density(spec,'polynomial',n0=n0,nedge=ne,exponent=6.0)
						profile = solver.N[spec][-1,:] / solver.N[spec][-1,0]
						plasma[i][j].set_density(spec,'interp',rho_vals=solver.rho_grid,n_vals=ntemp*profile)
						profile = solver.T[spec][-1,:] / solver.T[spec][-1,0]
						# plasma[i][j].set_temperature(spec,'polynomial',T0=t0,Tedge=te,exponent=1.0)
						plasma[i][j].set_temperature(spec,'interp',rho_vals=solver.rho_grid,T_vals=ttemp*profile)
			fren = tau_E[-1]/tau_ISS04[-1]
			popcon = POPCON(solver.B,solver.aminor,solver.Rmajor,solver.iota23,plasma, iss04_fact=fren, make_plot=False, popcon_title=f'fren={fren:.2f}')
			px = 1/plt.rcParams['figure.dpi']
			font = {'family' : 'Arial',
					'weight' : 'normal',
					'size'   : 24}
			matplotlib.rc('font', **font)
			fig,ax = plt.subplots(1,1,figsize=(1024*px,768*px))
			popcon.plot_popcon(ax=ax)
			ax.plot(solver.T['deuterium'][:,0]/1E3,solver.N['deuterium'][:,0]/1E20,color='black')
			ax.set_xlim([te_min/1E3,te_max/1E3])
			ax.set_ylim([ne_min/1E20,ne_max/1E20])
			plt.show()
			if (args.lsave): fig.savefig(f'popcon_{args.output_ext}.png', dpi=fig.dpi)
		if args.tslice_profs > 0:
			tdex = np.count_nonzero(solver.time<args.tslice_profs)
			px = 1/plt.rcParams['figure.dpi']
			font = {'family' : 'Arial',
					'weight' : 'normal',
					'size'   : 14}
			matplotlib.rc('font', **font)
			fig,ax = plt.subplots(4,2,figsize=(1024*px,1024*px))
			# Densities (main ion/electron)
			for spec in ['electrons','hydrogen','deuterium','tritium']:
				if spec in solver.list_of_species:
					ax[0,0].plot(solver.rho_grid,solver.N[spec][tdex,:]/1E19,label=spec,linewidth=2.0)
			ax[0,0].set_ylabel(r'Density $\times10^{19}$ [$m^{-3}$]')
			ax[0,0].set_xlim([0,1])
			ax[0,0].legend()
			# Temperatures
			for spec in ['electrons','hydrogen','deuterium','tritium']:
				if spec in solver.list_of_species:
					ax[0,1].plot(solver.rho_grid,solver.T[spec][tdex,:]/1000,label=spec,linewidth=2.0)
			ax[0,1].set_ylabel('Temperature [keV]')
			ax[0,1].set_xlim([0,1])
			ax[0,1].yaxis.set_label_position("right")
			ax[0,1].yaxis.tick_right()
			# Densities (main ion/electron)
			for spec in solver.list_of_species:
				if spec not in ['electrons','hydrogen','deuterium','tritium']:
					ax[1,0].plot(solver.rho_grid,solver.N[spec][tdex,:]/1E18,label=spec,linewidth=2.0)
			ax[1,0].set_ylabel(r'Density $\times10^{18}$ [$m^{-3}$]')
			ax[1,0].set_xlim([0,1])
			ax[1,0].legend()
			# Temperatures
			for spec in solver.list_of_species:
				if spec not in ['electrons','hydrogen','deuterium','tritium']:
					ax[1,1].plot(solver.rho_grid,solver.T[spec][tdex,:]/1000,label=spec,linewidth=2.0)
			ax[1,1].set_ylabel('Temperature [keV]')
			ax[1,1].set_xlim([0,1])
			ax[1,1].yaxis.set_label_position("right")
			ax[1,1].yaxis.tick_right()
			# Particle Sources (main ion/electron)
			for spec in ['electrons','hydrogen','deuterium','tritium']:
				if spec in solver.list_of_species:
					S = np.zeros((nt,nr))
					for source in solver.explicit_particle_sources[spec].keys():
						S += solver.explicit_particle_sources[spec][source]
					ax[2,0].plot(solver.rho_grid,S[tdex,:]/1E18,label=spec,linewidth=2.0)
			ax[2,0].set_ylabel(r'Source Rate $\times10^{18}$ [$m^{-3}/s$]')
			ax[2,0].set_xlim([0,1])
			ax[2,0].legend()
			# Heating Sources (main ion/electron)
			for spec in ['electrons','hydrogen','deuterium','tritium']:
				if spec in solver.list_of_species:
					S = np.zeros((nt,nr))
					for source in solver.explicit_energy_sources[spec].keys():
						S += solver.explicit_energy_sources[spec][source]
					ax[2,1].plot(solver.rho_grid,S[tdex,:]/1E3,label=spec,linewidth=2.0)
			ax[2,1].set_ylabel(r'Heating Rate [$kW/m^3$]')
			ax[2,1].set_xlim([0,1])
			ax[2,1].legend()
			ax[2,1].yaxis.set_label_position("right")
			ax[2,1].yaxis.tick_right()
			# Particle Sources (other)
			for spec in solver.list_of_species:
				if spec not in ['electrons','hydrogen','deuterium','tritium']:
					S = np.zeros((nt,nr))
					for source in solver.explicit_particle_sources[spec].keys():
						S += solver.explicit_particle_sources[spec][source]
					ax[3,0].plot(solver.rho_grid,S[tdex,:]/1E18,label=spec,linewidth=2.0)
			ax[3,0].set_ylabel(r'Source Rate $\times10^{18}$ [$m^{-3}/s$]')
			ax[3,0].set_xlabel(r'Radius (r/a)')
			ax[3,0].set_xlim([0,1])
			ax[3,0].legend()
			# Heating Sources (other)
			for spec in solver.list_of_species:
				if spec not in ['electrons','hydrogen','deuterium','tritium']:
					S = np.zeros((nt,nr))
					for source in solver.explicit_energy_sources[spec].keys():
						S += solver.explicit_energy_sources[spec][source]
					ax[3,1].plot(solver.rho_grid,S[tdex,:]/1E3,label=spec,linewidth=2.0)
			ax[3,1].set_ylabel(r'Heating Rate [$kW/m^3$]')
			ax[3,1].set_xlabel(r'Radius (r/a)')
			ax[3,1].set_xlim([0,1])
			ax[3,1].legend()
			ax[3,1].yaxis.set_label_position("right")
			ax[3,1].yaxis.tick_right()
			ax[0,0].set_title(f'Profiles at t={solver.time[tdex]}s')
			plt.show()
			if (args.lsave): fig.savefig(f'profs_{args.output_ext}_t{np.round(args.tslice_profs*1000)}ms.png', dpi=fig.dpi)


	sys.exit(0)
