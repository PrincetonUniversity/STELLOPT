#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Main routine
if __name__=="__main__":
	import sys
	from argparse import ArgumentParser
	import matplotlib.pyplot as pyplot
	from matplotlib.backends.backend_agg import FigureCanvasAgg
	from libstell.vmec import VMEC
	from libstell.nescoil import NESCOIL
	from libstell.fieldlines import FIELDLINES
	from libstell.plot3D import PLOT3D
	import numpy as np
	parser = ArgumentParser(description= 
		'''Provides a tool for accessing fieldlines data.''')
	parser.add_argument("-f", "--fieldlines", dest="fieldlines_ext",
		help="FIELDLINES file extension", default = None)
	parser.add_argument("-p", "--plot", dest="lplot", action='store_true',
		help="Plot the fieldlines file.", default = False)
	parser.add_argument("--plotphi", dest="plotphi",
		help="Poincare plot at a given index.", default = None, type=int)
	parser.add_argument("-i", "--info", dest="linfo", action='store_true',
		help="Print some info the the screen.", default = False)
	parser.add_argument("--plot3d", dest="k3d",
		help="Plot a fieldline in 3D.", default = None, type=int)
	parser.add_argument("--plotpoinc3d", dest="poinc3d",
		help="Plot a 3D Poincare plot.", default = None, type=int)
	parser.add_argument("--plotiota", dest="plot_iota", action='store_true',
		help="Plot the iota profile.", default = False)
	parser.add_argument("-v", "--vmec", dest="vmec_ext", 
		help="Add VMEC equilbrium to plot", default = None)
	parser.add_argument("--nescoil", dest="nescoil_file", 
		help="Add NESCOIL surfaces to the plot", default = None)
	parser.add_argument("--nskip", dest="nskip",
		help="Field line skipping parameter when generating Poincare cross sections (default: 1)", default = 1, type=int)
	parser.add_argument("--colormap", dest="colormap", 
		help="Colormap to use for plots (default: hot)", default = 'hot')
	parser.add_argument("--plotheat", dest="heatfactor",
		help="Plot the heatflux scaled to a total power in [W]", default = None, type=float)
	parser.add_argument("--plotindex", dest="i3d",
		help="Plot all fieldlines at a given Poincare index in 3D.", default = None, type=int)
	parser.add_argument("--plot_brz", dest="brz_index_phi",
		help="Plot the B-Field at an R/Z-plane, fixed phi (index)", default = None, type=int)
	parser.add_argument("--plot_brphi", dest="brphi_index_phi",
		help="Plot the B-Field at an R/phi-plane, fixed Z (index)", default = None, type=int)
	parser.add_argument("--plot_baxis", dest="lplot_baxis", action='store_true',
		help="Plots |B| along the magnetic axis (first field line).", default = False)
	parser.add_argument("--output_asc", dest="asc_phi",
		help="Output a given Poincare cross section at a given phi value [deg].", default = None, type=float)
	parser.add_argument("--output_asc_orbit", dest="asc_orbit",
		help="Output a given field line trajectory for a given fieldline.", default = None, type=int)
	parser.add_argument("--save", dest="lsave", action='store_true',
		help="Save the plots with ext names.", default = False)
	parser.add_argument("--background", dest="lbackground", action='store_true',
		help="Supress rendering window on plot.", default = False)
	args = parser.parse_args()
	field_data = FIELDLINES()
	px = 1/pyplot.rcParams['figure.dpi']
	if args.fieldlines_ext:
		field_data.read_fieldlines('fieldlines_'+args.fieldlines_ext+'.h5')
		if args.linfo:
			print(f'  NMARKERS:  {np.sum(field_data.nlines):d}')
			print(f'  WALL HITS: {np.sum(field_data.wall_strikes):d}')
		if type(args.plotphi) is not type(None):
			fig,ax = pyplot.subplots(1,1,figsize=(1024*px,768*px))
			if args.lbackground: canvas = FigureCanvasAgg(fig)
			phi0 = field_data.PHI_lines[0,args.plotphi]
			field_data.plot_poincare(phi0,args.nskip,ax=ax)
			if args.nescoil_file:
				nescout = NESCOIL()
				nescout.read_nescout(args.nescoil_file)
				theta = np.linspace([0],[2.0*np.pi],360)
				phi = np.array([[phi0*nescout.np]])
				nescout.computesurfaces(theta=theta,zeta=phi)
				ax1.plot(nescout.rp[0,:,0],nescout.zp[0,:,0],'r')
				ax1.plot(nescout.rc[0,:,0],nescout.zc[0,:,0],'b')
			if args.vmec_ext:
				vmec_wout = VMEC()
				vmec_wout.read_wout(args.vmec_ext)
				theta = np.ndarray((360,1))
				for j in range(360): theta[j]=2.0*np.pi*j/359.0
				phi = np.array([[phi0]])
				r = vmec_wout.cfunct(theta,phi,vmec_wout.rmnc,vmec_wout.xm,vmec_wout.xn)
				z = vmec_wout.sfunct(theta,phi,vmec_wout.zmns,vmec_wout.xm,vmec_wout.xn)
				j = vmec_wout.ns-1
				ax.plot(r[j,:,0],z[j,:,0],'r')
			if not args.lbackground:pyplot.show()
			if (args.lsave): fig.savefig(f'poincare_phi{args.plotphi:03d}_{args.fieldlines_ext}.png', dpi=fig.dpi)
		if args.lplot:
			fig,(ax1,ax2,ax3) = pyplot.subplots(1,3,sharey=True,figsize=(1024*px,512*px))
			if args.lbackground: canvas = FigureCanvasAgg(fig)
			pyplot.subplots_adjust(hspace=0.1,wspace=0.15)
			phi0 = 0
			field_data.plot_poincare(phi0,args.nskip,ax=ax1)
			phi1 = field_data.PHI_lines[0,int(np.round(field_data.npoinc/4))]
			field_data.plot_poincare(phi1,args.nskip,ax=ax2)
			phi2 = field_data.PHI_lines[0,int(np.round(field_data.npoinc/2))]
			field_data.plot_poincare(phi2,args.nskip,ax=ax3)
			if args.nescoil_file:
				nescout = NESCOIL()
				nescout.read_nescout(args.nescoil_file)
				theta = np.linspace([0],[2.0*np.pi],360)
				phi = np.array([[phi0*nescout.np],[phi1*nescout.np],[phi2*nescout.np]])
				nescout.computesurfaces(theta=theta,zeta=phi)
				ax1.plot(nescout.rp[0,:,0],nescout.zp[0,:,0],'r')
				ax2.plot(nescout.rp[0,:,1],nescout.zp[0,:,1],'r')
				ax3.plot(nescout.rp[0,:,2],nescout.zp[0,:,2],'r')
				ax1.plot(nescout.rc[0,:,0],nescout.zc[0,:,0],'b')
				ax2.plot(nescout.rc[0,:,1],nescout.zc[0,:,1],'b')
				ax3.plot(nescout.rc[0,:,2],nescout.zc[0,:,2],'b')
			if args.vmec_ext:
				vmec_wout = VMEC()
				vmec_wout.read_wout(args.vmec_ext)
				theta = np.ndarray((360,1))
				for j in range(360): theta[j]=2.0*np.pi*j/359.0
				phi = np.array([[phi0],[phi1],[phi2]])
				r = vmec_wout.cfunct(theta,phi,vmec_wout.rmnc,vmec_wout.xm,vmec_wout.xn)
				z = vmec_wout.sfunct(theta,phi,vmec_wout.zmns,vmec_wout.xm,vmec_wout.xn)
				j = vmec_wout.ns-1
				ax1.plot(r[j,:,0],z[j,:,0],'r')
				ax2.plot(r[j,:,1],z[j,:,1],'r')
				ax3.plot(r[j,:,2],z[j,:,2],'r')
			if not args.lbackground:pyplot.show()
			if (args.lsave): fig.savefig(f'poincare_{args.fieldlines_ext}.png', dpi=fig.dpi)
		if args.plot_iota:
			[r,iota,iota_err] = field_data.calc_iota()
			fig,ax = pyplot.subplots(1,1,figsize=(1024*px,768*px))
			ax.plot(r,iota,'ok',label='FIELDLINES')
			#ax.set_xlim([0,1])
			ax.set_ylim([0.75,1.25])
			ax.set_xlabel('Average Minor Radius [m]')
			ax.set_ylabel(r'$\iota$')
			if args.vmec_ext:
				vmec_wout = VMEC()
				vmec_wout.read_wout(args.vmec_ext)
				s = np.linspace(0,1.0,vmec_wout.ns)
				r = np.sqrt(s)*vmec_wout.aminor
				ax.plot(r,vmec_wout.iotaf,'r',linewidth=2.0,label='VMEC')
				pyplot.legend()
			ax.set_xlim([0,3.0])
			if not args.lbackground:pyplot.show()
			if (args.lsave): fig.savefig(f'iota_{args.fieldlines_ext}.png', dpi=fig.dpi)
		if args.poinc3d:
			plt3d = PLOT3D()
			field_data.plot_poincare3D(args.poinc3d,plot3D=plt3d,pointsize=0.1)
			plt3d.render(args.lbackground)
			if (args.lsave): plt3d.saveImage(f'poinc3d_{args.fieldlines_ext}.png')
		if args.k3d:
			plt3d = PLOT3D()
			field_data.plot_cloud(args.k3d,plot3D=plt3d,pointsize=0.1)
			plt3d.render(args.lbackground)
			if (args.lsave): plt3d.saveImage(f'poinc_cloud_{args.fieldlines_ext}.png')
		if args.i3d:
			plt3d = PLOT3D()
			field_data.plot_index3d(args.i3d,plot3D=plt3d,pointsize=0.1)
			plt3d.render(args.lbackground)
			if (args.lsave): plt3d.saveImage(f'poinc3d_{args.i3d:03d}_{args.fieldlines_ext}.png')
		if args.heatfactor:
			plt3d = PLOT3D()
			fact = args.heatfactor/field_data.nlines
			field_data.plot_heatflux(factor=args.heatfactor/field_data.nlines,colormap=args.colormap,plot3D=plt3d)
			plt3d.colorbar(title=rf'Q [W/$m^2$]')
			plt3d.render(args.lbackground)
			if (args.lsave): plt3d.saveImage(f'heatflux3d_{args.fieldlines_ext}.png')
		if type(args.brz_index_phi) is not type(None):
			fig,ax = pyplot.subplots(2,2,sharey=True,figsize=(1024*px,768*px))
			if args.lbackground: canvas = FigureCanvasAgg(fig)
			j = args.brz_index_phi
			x = np.squeeze(field_data.raxis)
			y = np.squeeze(field_data.zaxis)
			b = np.sqrt(field_data.B_R**2+field_data.B_Z**2+field_data.B_PHI**2)
			h0=ax[0,0].pcolormesh(x,y,np.squeeze(field_data.B_R[:,j,:]).T,cmap=args.colormap,shading='gouraud')
			ax[0,0].set_xlabel('R [m]'); ax[0,0].set_ylabel('Z [m]'); 
			h0.set_clim(vmin=-2.0,vmax=2.0); fig.colorbar(h0,label=r'$B_R$ [T]')
			h1=ax[0,1].pcolormesh(x,y,np.squeeze(field_data.B_Z[:,j,:]).T,cmap=args.colormap,shading='gouraud')
			ax[0,1].set_xlabel('R [m]'); ax[0,1].set_ylabel('Z [m]'); 
			h1.set_clim(vmin=-2.0,vmax=2.0); fig.colorbar(h1,label=r'$B_Z$ [T]')
			h2=ax[1,0].pcolormesh(x,y,np.squeeze(field_data.B_PHI[:,j,:]).T,cmap=args.colormap,shading='gouraud')
			ax[1,0].set_xlabel('R [m]'); ax[1,0].set_ylabel('Z [m]'); 
			h2.set_clim(vmin=-7.0,vmax=7.0); fig.colorbar(h2,label=r'$B_\phi$ [T]')
			h3=ax[1,1].pcolormesh(x,y,np.squeeze(b[:,j,:]).T,cmap=args.colormap,shading='gouraud')
			ax[1,1].set_xlabel('R [m]'); ax[1,1].set_ylabel('Z [m]'); 
			h3.set_clim(vmin=0.0,vmax=10.0); fig.colorbar(h3,label=r'$|B|$ [T]')
			if not args.lbackground:pyplot.show()
			if (args.lsave): fig.savefig(f'brz_{args.brz_index_phi:0.3d}_{args.fieldlines_ext}.png', dpi=fig.dpi)
		if type(args.brphi_index_phi) is not type(None):
			fig,ax = pyplot.subplots(2,2,sharey=True,figsize=(1024*px,768*px))
			if args.lbackground: canvas = FigureCanvasAgg(fig)
			j = args.brphi_index_phi
			x = np.squeeze(field_data.raxis)
			y = np.squeeze(field_data.phiaxis)
			b = np.sqrt(field_data.B_R**2+field_data.B_Z**2+field_data.B_PHI**2)
			h0=ax[0,0].pcolormesh(x,y,np.squeeze(field_data.B_R[:,:,j]).T,cmap=args.colormap,shading='gouraud')
			ax[0,0].set_xlabel('R [m]'); ax[0,0].set_ylabel(r'$\phi$ [rad]'); 
			h0.set_clim(vmin=-1.0,vmax=1.0); fig.colorbar(h0,label=r'$B_R$ [T]')
			h1=ax[0,1].pcolormesh(x,y,np.squeeze(field_data.B_Z[:,:,j]).T,cmap=args.colormap,shading='gouraud')
			ax[0,1].set_xlabel('R [m]'); ax[0,1].set_ylabel(r'$\phi$ [rad]'); 
			h1.set_clim(vmin=-1.0,vmax=1.0); fig.colorbar(h1,label=r'$B_Z$ [T]')
			h2=ax[1,0].pcolormesh(x,y,np.squeeze(field_data.B_PHI[:,:,j]).T,cmap=args.colormap,shading='gouraud')
			ax[1,0].set_xlabel('R [m]'); ax[1,0].set_ylabel(r'$\phi$ [rad]'); 
			h2.set_clim(vmin=-7.0,vmax=7.0); fig.colorbar(h2,label=r'$B_\phi$ [T]')
			h3=ax[1,1].pcolormesh(x,y,np.squeeze(b[:,:,j]).T,cmap=args.colormap,shading='gouraud')
			ax[1,1].set_xlabel('R [m]'); ax[1,1].set_ylabel(r'$\phi$ [rad]'); 
			h3.set_clim(vmin=0.0,vmax=10.0); fig.colorbar(h3,label=r'$|B|$ [T]')
			if not args.lbackground:pyplot.show()
			if (args.lsave): fig.savefig(f'brphi_{args.brphi_index_phi:0.3d}_{args.fieldlines_ext}.png', dpi=fig.dpi)
		if args.lplot_baxis:
			fig,ax = pyplot.subplots(1,1,figsize=(1024*px,768*px))
			if args.lbackground: canvas = FigureCanvasAgg(fig)
			phi = np.rad2deg(field_data.PHI_lines[0,:-2])
			high = np.rad2deg(field_data.phiaxis[-1]/2.0)
			low  = -high
			phi = ((phi - low) % (high-low)) + low
			ax.plot(phi,field_data.B_lines[0,:-2],'k.',linewidth=4)
			ax.set_title('|B| along magnetic axis')
			ax.set_xlabel(r'Toroidal Angle $\phi$ [$^o$]')
			ax.set_ylabel('|B| [T]')
			if not args.lbackground:pyplot.show()
			if (args.lsave): fig.savefig(f'baxis_{args.fieldlines_ext}.png', dpi=fig.dpi)
		if type(args.asc_phi) is not type(None):
			field_data.write_asc([np.deg2rad(args.asc_phi)],nskip=args.nskip,filename=f'poincare_{args.fieldlines_ext}_phi_{int(args.asc_phi):03d}.asc')
		if type(args.asc_orbit) is not type(None):
			field_data.write_orbit_asc(args.asc_orbit,filename=f'poincare_{args.fieldlines_ext}_k_{int(args.asc_orbit):03d}.asc')
	sys.exit(0)