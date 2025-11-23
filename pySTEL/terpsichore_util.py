#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Main routine
if __name__=="__main__":
	import sys
	from argparse import ArgumentParser
	import numpy as np
	import matplotlib.pyplot as pyplot
	from libstell.vmec import VMEC
	from libstell.terpsichore import TERPSICHORE
	parser = ArgumentParser(description= 
		'''Provides a simple tool for creating and plotting TERPSICHORE files.''')
	parser.add_argument("--vmec", dest="vmec_ext",
		help="VMEC wout file", default = None)
	parser.add_argument("--terpfile", dest="terp_ext",
		help="TERPSICHORE File", default = None)
	parser.add_argument("--input", dest="linput", action='store_true',
		help="Create input files for TEPRSICHORE.", default = False)
	parser.add_argument("--plot_input", dest="lplot_in", action='store_true',
		help="Plot the wall data in fort.17.", default = False)
	parser.add_argument("--plot_wall", dest="lplot_wall", action='store_true',
		help="Plot the wall data in fort.19.", default = False)
	parser.add_argument("--plot_22", dest="lplot_22", action='store_true',
		help="Plot the wall data in fort.22.", default = False)
	parser.add_argument("--plot_23", dest="lplot_23", action='store_true',
		help="Plot the wall data in fort.23.", default = False)
	terp_data = TERPSICHORE()
	vmec_data = VMEC()
	args = parser.parse_args()
	if args.linput and args.vmec_ext:
		try:
			vmec_data.read_wout(args.vmec_ext)
		except:
			print(f'Could not open VMEC wout file: {args.vmec_ext}')
			sys.exit(-1)
		terp_data.write_eq_input(vmec_data)
		nmodes = int(np.floor(vmec_data.nfp/2.0))
		print('Values for tpr_modules_ap.f')
		for n in range(nmodes+1):
			terp_data.create_input(vmec_data,n)
	if args.lplot_in:
		if args.terp_ext:
			file = args.terp_ext
		else:
			file = 'fort.17'
		terp_data.read_terpsichore_17(file)
		theta = np.linspace([0.0],[2*np.pi],128)
		zeta  = np.linspace([0.0],[  np.pi]/terp_data.nfp_in,3)
		r = terp_data.cfunct(theta,zeta,terp_data.rmnc_in,terp_data.xm_in,terp_data.xn_in)
		z = terp_data.sfunct(theta,zeta,terp_data.zmns_in,terp_data.xm_in,terp_data.xn_in)
		px = 1/pyplot.rcParams['figure.dpi']
		fig=pyplot.figure(figsize=(1024*px,768*px))
		ax=fig.add_subplot(131)
		ax.plot(np.squeeze(r[:,:,0]).T,np.squeeze(z[:,:,0]).T,'k')
		ax.set_xlabel('R [m]'); ax.set_ylabel('Z [m]')
		ax.set_aspect('equal', adjustable='box')
		ax=fig.add_subplot(132)
		ax.plot(np.squeeze(r[:,:,1]).T,np.squeeze(z[:,:,1]).T,'k')
		ax.set_xlabel('R [m]'); ax.set_ylabel('Z [m]')
		ax.set_aspect('equal', adjustable='box')
		ax=fig.add_subplot(133)
		ax.plot(np.squeeze(r[:,:,2]).T,np.squeeze(z[:,:,2]).T,'k')
		ax.set_xlabel('R [m]'); ax.set_ylabel('Z [m]')
		ax.set_aspect('equal', adjustable='box')
		pyplot.show()
	if args.lplot_wall:
		if args.terp_ext:
			file = args.terp_ext
		else:
			file = 'fort.19'
		terp_data.read_terpsichore_19(file)
		px = 1/pyplot.rcParams['figure.dpi']
		fig=pyplot.figure(figsize=(1024*px,768*px))
		ax=fig.add_subplot(131)
		ax.plot(np.squeeze(terp_data.Rpvi[:,0]),np.squeeze(terp_data.Zpvi[:,0]),'r')
		ax.plot(np.squeeze(terp_data.Rwall[:,0]),np.squeeze(terp_data.Zwall[:,0]),'k')
		ax.set_xlabel('R [m]'); ax.set_ylabel('Z [m]')
		ax.set_aspect('equal', adjustable='box')
		ax=fig.add_subplot(132)
		j = round(terp_data.nk/4)
		ax.plot(np.squeeze(terp_data.Rpvi[:,j]),np.squeeze(terp_data.Zpvi[:,j]),'r')
		ax.plot(np.squeeze(terp_data.Rwall[:,j]),np.squeeze(terp_data.Zwall[:,j]),'k')
		ax.set_xlabel('R [m]'); ax.set_ylabel('Z [m]')
		ax.set_aspect('equal', adjustable='box')
		ax=fig.add_subplot(133)
		j = round(terp_data.nk/2)
		ax.plot(np.squeeze(terp_data.Rpvi[:,j]),np.squeeze(terp_data.Zpvi[:,j]),'r')
		ax.plot(np.squeeze(terp_data.Rwall[:,j]),np.squeeze(terp_data.Zwall[:,j]),'k')
		ax.set_xlabel('R [m]'); ax.set_ylabel('Z [m]')
		ax.set_aspect('equal', adjustable='box')
		pyplot.show()
	if args.lplot_22:
		if args.terp_ext:
			file = arg.terp_ext
		else:
			file = 'fort.22'
		terp_data.read_terpsichore_22(file)
		s_terp = np.linspace(0.0,1.0,terp_data.ni)
		px = 1/pyplot.rcParams['figure.dpi']
		fig=pyplot.figure(figsize=(1024*px,768*px))
		ax=fig.add_subplot(231)
		ax.plot(s_terp,terp_data.pth/(np.pi*4E-7),'k')
		ax.set_xlabel('Radial Grid')
		ax.set_ylabel('Pressure [Pa]')
		ax.set_title('Thermal Pressure')
		ax.set_xlim([0,1])
		ax.text(0.25,0.9*max(terp_data.pth/(np.pi*4E-7)),f'WP/WK = {terp_data.gamma}')
		ax=fig.add_subplot(232)
		ax.plot(s_terp,terp_data.aiota,'k')
		ax.set_xlabel('Radial Grid')
		ax.set_ylabel('Iota')
		ax.set_title('Rotational Transform')
		ax.set_xlim([0,1])
		ax=fig.add_subplot(233)
		ax.plot(s_terp,terp_data.wpsi,'k')
		ax.set_xlabel('Radial Grid')
		ax.set_ylabel(r'$\delta W$')
		ax.set_title('Perturbed Energy')
		ax.set_xlim([0,1])
		ax=fig.add_subplot(234)
		ax.plot(s_terp,terp_data.am.T)
		ax.set_xlabel('Radial Grid')
		ax.set_ylabel('AM')
		ax.set_title('Perturbed Pressure Modes')
		ax.set_xlim([0,1])
		ax=fig.add_subplot(235)
		ax.plot(s_terp,terp_data.pvp.T)
		ax.set_xlabel('Radial Grid')
		ax.set_ylabel(r'$\eta$ (Binormal)')
		ax.set_title(r'Perturbed $\eta$ Modes')
		ax.set_xlim([0,1])
		ax=fig.add_subplot(236)
		ax.plot(s_terp,terp_data.pvpi.T)
		ax.set_xlabel('Radial Grid')
		ax.set_ylabel(r'$\mu$ (Parallel)')
		ax.set_title(r'Perturbed $\mu$ Modes')
		ax.set_xlim([0,1])
		pyplot.subplots_adjust(wspace=0.50,hspace=0.4)
		pyplot.show()
	if args.lplot_23:
		if args.terp_ext:
			file = arg.terp_ext
		else:
			file = 'fort.23'
		terp_data.read_terpsichore_23(file)
		# Plot all modes
		px = 1/pyplot.rcParams['figure.dpi']
		fig=pyplot.figure(figsize=(1024*px,768*px))
		ax=fig.add_subplot(121)
		ax.plot(terp_data.s[0:],terp_data.xi)
		ax.set_xlabel('Norm. Toroidal Flux (s)'); ax.set_ylabel(r'$\xi$ Modes')
		ax.set_title('Perturbed Radial Modes')
		ax.text(0.95,0.98,rf'$W_p$={terp_data.wp:12.3E}', horizontalalignment='right',\
				verticalalignment='center', transform=ax.transAxes)
		ax.text(0.95,0.93,rf'$W_k$={terp_data.wk:12.3E}', horizontalalignment='right',\
				verticalalignment='center', transform=ax.transAxes)
		ax.text(0.95,0.88,rf'$W_p/W_k$={terp_data.wp/terp_data.wk:12.3E}', horizontalalignment='right',\
				verticalalignment='center', transform=ax.transAxes)
		if (terp_data.wp/terp_data.wk < 0):
			ax.text(0.95,0.83,rf'Unstable', horizontalalignment='right',\
					verticalalignment='center', transform=ax.transAxes)
		else:
			ax.text(0.95,0.83,rf'Stable', horizontalalignment='right',\
					verticalalignment='center', transform=ax.transAxes)
		ax=fig.add_subplot(122)
		ax.plot(terp_data.s[1:],terp_data.eta)
		ax.set_xlabel('Norm. Toroidal Flux (s)'); ax.set_ylabel(rf'$\eta$ Modes')
		ax.set_title('Perturbed Binormal Modes')
		pyplot.show()
		# Plot largest 5 modes
		ntop = 5
		eta_max = np.amax(terp_data.eta,axis=0)
		xi_max  = np.amax(terp_data.xi,axis=0)
		eta_ind = np.argpartition(eta_max, -ntop)[-ntop:]
		xi_ind  = np.argpartition(xi_max,  -ntop)[-ntop:]
		eta_ind = eta_ind[np.argsort(eta_max[eta_ind])]
		xi_ind  = xi_ind[np.argsort(xi_max[xi_ind])]
		px = 1/pyplot.rcParams['figure.dpi']
		fig=pyplot.figure(figsize=(1024*px,768*px))
		ax=fig.add_subplot(121)
		for i in xi_ind[::-1]:
			m = int(terp_data.mb[i,0])
			n = int(terp_data.nb[i,0]/terp_data.nper)
			ax.plot(terp_data.s[0:],terp_data.xi[:,i],label=rf'm/n={m:2d}/{n:2d}')
		ax.set_xlabel('Norm. Toroidal Flux (s)'); ax.set_ylabel(rf'$\xi$ Modes')
		ax.set_title('Perturbed Radial Modes (5 Largest)')
		ax.legend()
		ax=fig.add_subplot(122)
		for i in eta_ind[::-1]:
			m = int(terp_data.mb[i,0])
			n = int(terp_data.nb[i,0]/terp_data.nper)
			ax.plot(terp_data.s[1:],terp_data.eta[:,i],label=rf'm/n={m:2d}/{n:2d}')
		ax.set_xlabel('Norm. Toroidal Flux (s)'); ax.set_ylabel(rf'$\eta$ Modes')
		ax.set_title('Perturbed Binormal Modes (5 Largest)')
		ax.legend()
		pyplot.show()

	sys.exit(0)