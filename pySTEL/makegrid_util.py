#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Main routine
if __name__=="__main__":
	import sys
	import os.path
	from argparse import ArgumentParser
	import matplotlib.pyplot as pyplot
	from matplotlib.backends.backend_agg import FigureCanvasAgg
	from libstell.vmec import VMEC, VMEC_INDATA
	from libstell.makegrid import MAKEGRID
	import numpy as np
	parser = ArgumentParser(description= 
		'''Provides a tool for accessing makegrid data.''')
	parser.add_argument("-m", "--makegrid", dest="mgrid_ext",
		help="Makegrid (mgrid) file extension", default = None)
	parser.add_argument("-v", "--vmec", dest="vmec_ext", 
		help="Add VMEC equilbrium to plot", default = None)
	parser.add_argument("--plot_brz", dest="brz_index_phi",
		help="Plot the B-Field at an R/Z-plane, fixed phi (index)", default = None, type=int)
	parser.add_argument("--plot_brphi", dest="brphi_index_phi",
		help="Plot the B-Field at an R/phi-plane, fixed Z (index)", default = None, type=int)
	parser.add_argument("--colormap", dest="colormap", 
		help="Colormap to use for plots (default: hot)", default = 'hot')
	parser.add_argument("--save", dest="lsave", action='store_true',
		help="Save the plots with ext names.", default = False)
	parser.add_argument("--background", dest="lbackground", action='store_true',
		help="Supress rendering window on plot.", default = False)
	args = parser.parse_args()
	mgrid_data = MAKEGRID()
	vmec_wout = VMEC()
	vmec_input = VMEC_INDATA()
	px = 1/pyplot.rcParams['figure.dpi']
	# Try to find VMEC file
	if type(args.vmec_ext) == type(None):
		print('  VMEC and MGRID file required.')
		sys.exit(-1)
	if os.path.isfile('input.'+args.vmec_ext):
		vmec_input.read_indata('input.'+args.vmec_ext)
		nv  = vmec_input.nzeta
		nfp = vmec_input.nfp
		extcur = vmec_input.extcur
	else:
		print(f'Could not find input file: input.{args.vmec_ext}')
		sys.exit(-1)
	# Now read file
	if args.mgrid_ext:
		mgrid_data.read_mgrid('mgrid_'+args.mgrid_ext+'.nc',extcur,nv,nfp)
		if type(args.brz_index_phi) is not type(None):
			fig,ax = pyplot.subplots(2,2,sharey=True,figsize=(1024*px,768*px))
			if args.lbackground: canvas = FigureCanvasAgg(fig)
			j = args.brz_index_phi
			x = np.linspace(mgrid_data.rminb,mgrid_data.rmaxb,mgrid_data.nr0b)
			y = np.linspace(mgrid_data.zminb,mgrid_data.zmaxb,mgrid_data.nz0b)
			b = np.sqrt(mgrid_data.brvac**2+mgrid_data.bzvac**2+mgrid_data.bpvac**2)
			h0=ax[0,0].pcolormesh(x,y,np.squeeze(mgrid_data.brvac[:,:,j]).T,cmap=args.colormap,shading='gouraud')
			ax[0,0].set_xlabel('R [m]'); ax[0,0].set_ylabel('Z [m]'); 
			h0.set_clim(vmin=-2.0,vmax=2.0); fig.colorbar(h0,label=r'$B_R$ [T]')
			h1=ax[0,1].pcolormesh(x,y,np.squeeze(mgrid_data.bzvac[:,:,j]).T,cmap=args.colormap,shading='gouraud')
			ax[0,1].set_xlabel('R [m]'); ax[0,1].set_ylabel('Z [m]'); 
			h1.set_clim(vmin=-2.0,vmax=2.0); fig.colorbar(h1,label=r'$B_Z$ [T]')
			h2=ax[1,0].pcolormesh(x,y,np.squeeze(mgrid_data.bpvac[:,:,j]).T,cmap=args.colormap,shading='gouraud')
			ax[1,0].set_xlabel('R [m]'); ax[1,0].set_ylabel('Z [m]'); 
			h2.set_clim(vmin=-7.0,vmax=7.0); fig.colorbar(h2,label=r'$B_\phi$ [T]')
			h3=ax[1,1].pcolormesh(x,y,np.squeeze(b[:,:,j]).T,cmap=args.colormap,shading='gouraud')
			ax[1,1].set_xlabel('R [m]'); ax[1,1].set_ylabel('Z [m]'); 
			h3.set_clim(vmin=0.0,vmax=10.0); fig.colorbar(h3,label=r'$|B|$ [T]')
			if not args.lbackground:pyplot.show()
			if (args.lsave): fig.savefig(f'brz_{args.brz_index_phi:0.3d}_{args.fieldlines_ext}.png', dpi=fig.dpi)
		if type(args.brphi_index_phi) is not type(None):
			fig,ax = pyplot.subplots(2,2,sharey=True,figsize=(1024*px,768*px))
			if args.lbackground: canvas = FigureCanvasAgg(fig)
			j = args.brphi_index_phi
			x = np.linspace(mgrid_data.rminb,mgrid_data.rmaxb,mgrid_data.nr0b)
			y = np.linspace(0,2.0*np.pi/nfp,nv,endpoint=False)
			b = np.sqrt(field_data.B_R**2+field_data.B_Z**2+field_data.B_PHI**2)
			h0=ax[0,0].pcolormesh(x,y,np.squeeze(mgrid_data.brvac[:,j,:]).T,cmap=args.colormap,shading='gouraud')
			ax[0,0].set_xlabel('R [m]'); ax[0,0].set_ylabel(r'$\phi$ [rad]'); 
			h0.set_clim(vmin=-1.0,vmax=1.0); fig.colorbar(h0,label=r'$B_R$ [T]')
			h1=ax[0,1].pcolormesh(x,y,np.squeeze(mgrid_data.bzvac[:,j,:]).T,cmap=args.colormap,shading='gouraud')
			ax[0,1].set_xlabel('R [m]'); ax[0,1].set_ylabel(r'$\phi$ [rad]'); 
			h1.set_clim(vmin=-1.0,vmax=1.0); fig.colorbar(h1,label=r'$B_Z$ [T]')
			h2=ax[1,0].pcolormesh(x,y,np.squeeze(mgrid_data.bpvac[:,j,:]).T,cmap=args.colormap,shading='gouraud')
			ax[1,0].set_xlabel('R [m]'); ax[1,0].set_ylabel(r'$\phi$ [rad]'); 
			h2.set_clim(vmin=-7.0,vmax=7.0); fig.colorbar(h2,label=r'$B_\phi$ [T]')
			h3=ax[1,1].pcolormesh(x,y,np.squeeze(b[:,j,:]).T,cmap=args.colormap,shading='gouraud')
			ax[1,1].set_xlabel('R [m]'); ax[1,1].set_ylabel(r'$\phi$ [rad]'); 
			h3.set_clim(vmin=0.0,vmax=10.0); fig.colorbar(h3,label=r'$|B|$ [T]')
			if not args.lbackground:pyplot.show()
			if (args.lsave): fig.savefig(f'brphi_{args.brphi_index_phi:0.3d}_{args.fieldlines_ext}.png', dpi=fig.dpi)
	sys.exit(0)