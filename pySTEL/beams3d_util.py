#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Main routine
if __name__=="__main__":
	import sys
	from argparse import ArgumentParser
	import matplotlib.pyplot as pyplot
	from libstell.vmec import VMEC
	from libstell.beams3d import BEAMS3D
	from libstell.plot3D import PLOT3D
	from libstell.wall import WALL
	from datetime import datetime
	import numpy as np
	parser = ArgumentParser(description= 
		'''Provides a tool for accessing beamsed data.''')
	parser.add_argument("-b", "--beams3d", dest="beams3d_ext",
		help="BEAMS3D file extension", default = None)
	parser.add_argument("-v", "--vmec", dest="vmec_ext", 
		help="Add VMEC equilbrium to plot", default = None)
	parser.add_argument("--plotloss", dest="lplotloss", action='store_true',
		help="Plot the loss rate.", default = False)
	parser.add_argument("--plotheat", dest="lplotheat", action='store_true',
		help="Plot the wall heatflux.", default = False)
	parser.add_argument("--plotorbits", dest="lplotorbits", action='store_true',
		help="Plot the 3D orbits.", default = False)
	parser.add_argument("--plotshine", dest="lplotshine", action='store_true',
		help="Plot the wall shinethrough.", default = False)
	parser.add_argument("--plotinjection", dest="lplotinjection", action='store_true',
		help="Plot the NBI injection.", default = False)
	parser.add_argument("--plottransport", dest="lplottrans", action='store_true',
		help="Plot the transport quantities.", default = False)
	parser.add_argument("--plotindex", dest="i3d",
		help="Plot all markers at a given orbit index in 3D.", default = None, type=int)
	parser.add_argument('--beams', nargs='+', dest="beams",
		help="List of beams to include.", default = None, type=int)
	parser.add_argument("--plotwall", dest="lplotwall", action='store_true',
		help="Plot the wall data if orbit plots.", default = False)
	parser.add_argument("--plot_surz", dest="su_index_phi",
		help="Plot the S and U field at an R/Z-plane, fixed phi (index)", default = None, type=int)
	parser.add_argument("--plot_brz", dest="brz_index_phi",
		help="Plot the B-Field at an R/Z-plane, fixed phi (index)", default = None, type=int)
	parser.add_argument("--plot_brphi", dest="brphi_index_phi",
		help="Plot the B-Field at an R/phi-plane, fixed Z (index)", default = None, type=int)
	parser.add_argument("--save", dest="lsave", action='store_true',
		help="Save the plots with ext names.", default = False)
	parser.add_argument("--background", dest="lbackground", action='store_true',
		help="Supress rendering window on plot.", default = False)
	args = parser.parse_args()
	beam_data = BEAMS3D()
	px = 1/pyplot.rcParams['figure.dpi']
	wall = None
	if args.vmec_ext:
		vmec_data=VMEC()
		vmec_data.read_wout('wout_'+args.vmec_ext+'.nc')
	if args.beams3d_ext:
		beam_data.read_beams3d('beams3d_'+args.beams3d_ext+'.h5')
		if args.lplotwall:
			wall = WALL()
			wall.name = 'BEAMS3D_WALL'
			wall.date = datetime.today().strftime('%Y-%m-%d')
			wall.nfaces = beam_data.nface
			wall.nvertex = beam_data.nvertex
			wall.vertex = beam_data.wall_vertex
			wall.faces = beam_data.wall_faces
			wall.laccel = False
		if args.lplotloss:
			time,nlost=beam_data.calcLoss()
			px = 1/pyplot.rcParams['figure.dpi']
			fig=pyplot.figure(figsize=(1024*px,768*px))
			ax=fig.add_subplot(111)
			for i in range(beam_data.nbeams):
				dex = beam_data.Beam == (i+1)
				npart = np.sum(beam_data.Weight[dex])
				ax.plot(time,100.*nlost[i,:]/npart,label=rf'Beam {int(i+1)}')
			ax.set_xscale('log')
			ax.set_xlabel('Time [s]')
			ax.set_ylabel('Particles Lost [%]')
			ax.set_title('Particle Loss Evolution')
			ax.legend()
			if not args.lbackground:pyplot.show()
			if (args.lsave): fig.savefig(f'loss_{args.beams3d_ext}.png', dpi=fig.dpi)
		if args.lplotheat:
			plt3d = PLOT3D()
			beam_data.plot_heatflux(beams=args.beams,load_type='heatflux',colormap='hot',plot3D=plt3d)
			plt3d.setClim(0,10E6)
			plt3d.colorbar(title=rf'Q [W/$m^2$]')
			plt3d.render(args.lbackground)
			if (args.lsave): plt3d.saveImage(f'heatflux_{args.beams3d_ext}.png')
		if args.lplotorbits:
			plt3d = PLOT3D()
			beam_data.plotorbit(plot3D=plt3d)
			if args.lplotwall:
				nfp = int(2*np.pi/beam_data.phiaxis[-1])
				zeta = np.linspace(0,np.pi,4*nfp+1)
				for z in zeta:
					wall.plot_wall_3D_RZ(z,plot3D=plt3d)
				wall.plot_wall_3D_XY(0.0,plot3D=plt3d)
			plt3d.render(args.lbackground)
			if (args.lsave): plt3d.saveImage(f'orbit_{args.beams3d_ext}.png')
		if args.lplotinjection:
			plt3d = PLOT3D()
			# Define markers and set colors
			mark = np.flatnonzero(beam_data.end_state == 0)
			beam_data.plotorbit(plot3D=plt3d,markers=mark,color='green')
			mark = np.flatnonzero(beam_data.end_state == 3)
			beam_data.plotorbit(plot3D=plt3d,markers=mark,color='red')
			mark = np.flatnonzero(beam_data.end_state == 4)
			beam_data.plotorbit(plot3D=plt3d,markers=mark,color='red')
			if args.lplotwall:
				nfp = int(2*np.pi/beam_data.phiaxis[-1])
				zeta = np.linspace(0,np.pi,4*nfp+1)
				for z in zeta:
					wall.plot_wall_3D_RZ(z,plot3D=plt3d)
				wall.plot_wall_3D_XY(0.0,plot3D=plt3d)
			plt3d.render(args.lbackground)
			if (args.lsave): plt3d.saveImage(f'injection_{args.beams3d_ext}.png')
		if args.lplotshine:
			plt3d = PLOT3D()
			beam_data.plot_heatflux(beams=args.beams,load_type='shine',colormap='hot',plot3D=plt3d)
			plt3d.setClim(0,10E6)
			plt3d.colorbar(title=rf'Q [W/$m^2$]')
			plt3d.render(args.lbackground)
			if (args.lsave): plt3d.saveImage(f'shine_{args.beams3d_ext}.png')
		if args.lplottrans:
			px = 1/pyplot.rcParams['figure.dpi']
			fig=pyplot.figure(figsize=(1024*px,768*px))
			ax=fig.add_subplot(221)
			ax2 = ax.twinx()
			pyplot.subplots_adjust(hspace=0.4,wspace=0.3)
			s,ne,ni,te,ti,zeff=beam_data.calcProfiles()
			ax.plot(np.sqrt(s),ne/1.0E19,'k')
			for i in range(ni.shape[0]):
				ax.plot(np.sqrt(s),ni[i,:]/1.0E19)
			ax2.plot(np.sqrt(s),te/1.0E3,'r')
			ax2.plot(np.sqrt(s),ti/1.0E3,'r--')
			ax.set_xlim(0,1)
			ax.set_xlabel('r/a')
			ax.set_ylabel(r'Density $x10^{19}$ [$m^{-3}$] ',color='k')
			ax2.set_ylabel(r"Temperature [keV]",color='r')
			ax=fig.add_subplot(222)
			ax2 = ax.twinx()
			aminor = beam_data.calcAminor()
			s, plasma_vol, plasma_dvolds = beam_data.calcVolume(ns=127)
			s, V, dVds=beam_data.calcEr()
			ax.plot(np.sqrt(s),V/1.0E3,'k')
			ax.set_xlabel('r/a')
			ax.set_ylabel(r'Potential [kV]',color='k')
			#ax2.set_ylabel(r'E_r (dV/ds) [kV/m^2]',color='r')
			ax=fig.add_subplot(223)
			s, births = beam_data.calcDepo(ns=127,beams=args.beams)
			ax.plot(np.sqrt(s),np.sum(births,axis=0)/1E19,'k')
			ax.set_xlabel('r/a')
			ax.set_ylabel(r'Birth Rate x10^{19} [$part/m^{-3}s$]')
			if not args.lbackground:pyplot.show()
			if (args.lsave): fig.savefig(f'transport_{args.beams3d_ext}.png', dpi=fig.dpi)
		if args.i3d:
			plt3d = PLOT3D()
			beam_data.plot_index3d(args.i3d,plot3D=plt3d,pointsize=0.1)
			if args.lplotwall:
				nfp = int(2*np.pi/beam_data.phiaxis[-1])
				zeta = np.linspace(0,np.pi,4*nfp+1)
				for z in zeta:
					wall.plot_wall_3D_RZ(z,plot3D=plt3d)
				wall.plot_wall_3D_XY(0.0,plot3D=plt3d)
			plt3d.render(args.lbackground)
			if (args.lsave): plt3d.saveImage(f'timestep_{args.i3d:03d}_{args.beams3d_ext}.png')
		if type(args.brz_index_phi) is not type(None):
			fig,ax = pyplot.subplots(2,2,sharey=True,figsize=(1024*px,768*px))
			j = args.brz_index_phi
			x = np.squeeze(beam_data.raxis)
			y = np.squeeze(beam_data.zaxis)
			b = np.sqrt(beam_data.B_R**2+beam_data.B_Z**2+beam_data.B_PHI**2)
			h0=ax[0,0].pcolormesh(x,y,np.squeeze(beam_data.B_R[:,j,:]).T,cmap='jet',shading='gouraud')
			ax[0,0].set_xlabel('R [m]'); ax[0,0].set_ylabel('Z [m]'); 
			h0.set_clim(vmin=-2.0,vmax=2.0); fig.colorbar(h0,label=r'$B_R$ [T]')
			h1=ax[0,1].pcolormesh(x,y,np.squeeze(beam_data.B_Z[:,j,:]).T,cmap='jet',shading='gouraud')
			ax[0,1].set_xlabel('R [m]'); ax[0,1].set_ylabel('Z [m]'); 
			h1.set_clim(vmin=-2.0,vmax=2.0); fig.colorbar(h1,label=r'$B_Z$ [T]')
			h2=ax[1,0].pcolormesh(x,y,np.squeeze(beam_data.B_PHI[:,j,:]).T,cmap='jet',shading='gouraud')
			ax[1,0].set_xlabel('R [m]'); ax[1,0].set_ylabel('Z [m]'); 
			h2.set_clim(vmin=-7.0,vmax=7.0); fig.colorbar(h2,label=r'$B_\phi$ [T]')
			h3=ax[1,1].pcolormesh(x,y,np.squeeze(b[:,j,:]).T,cmap='jet',shading='gouraud')
			ax[1,1].set_xlabel('R [m]'); ax[1,1].set_ylabel('Z [m]'); 
			h3.set_clim(vmin=0.0,vmax=10.0); fig.colorbar(h3,label=r'$|B|$ [T]')
			if not args.lbackground:pyplot.show()
			if (args.lsave): fig.savefig(f'brz_{args.brz_index_phi:03d}_{args.beams3d_ext}.png', dpi=fig.dpi)
		if type(args.su_index_phi) is not type(None):
			fig,ax = pyplot.subplots(1,2,sharey=True,figsize=(1024*px,768*px))
			j = args.su_index_phi
			x = np.squeeze(beam_data.raxis)
			y = np.squeeze(beam_data.zaxis)
			h0=ax[0].pcolormesh(x,y,np.squeeze(beam_data.S_ARR[:,j,:]).T,cmap='jet',shading='gouraud')
			ax[0].set_xlabel('R [m]'); ax[0].set_ylabel('Z [m]'); 
			h0.set_clim(vmin=0.0,vmax=4.0); fig.colorbar(h0,label=r'$Norm. Tor. Flux$ (s) [arb]')
			h1=ax[1].pcolormesh(x,y,np.squeeze(beam_data.U_ARR[:,j,:]).T,cmap='jet',shading='gouraud')
			ax[1].set_xlabel('R [m]'); ax[1].set_ylabel('Z [m]'); 
			h1.set_clim(vmin=0,vmax=2.*np.pi); fig.colorbar(h1,label=r'$Pol. Angle$ [rad]')
			if not args.lbackground:pyplot.show()
			if (args.lsave): fig.savefig(f'surz_{args.su_index_phi:03d}_{args.beams3d_ext}.png', dpi=fig.dpi)
		if type(args.brphi_index_phi) is not type(None):
			fig,ax = pyplot.subplots(2,2,sharey=True,figsize=(1024*px,768*px))
			j = args.brphi_index_phi
			x = np.squeeze(beam_data.raxis)
			y = np.squeeze(beam_data.phiaxis)
			b = np.sqrt(beam_data.B_R**2+beam_data.B_Z**2+beam_data.B_PHI**2)
			h0=ax[0,0].pcolormesh(x,y,np.squeeze(beam_data.B_R[:,:,j]).T,cmap='jet',shading='gouraud')
			ax[0,0].set_xlabel('R [m]'); ax[0,0].set_ylabel(r'$\phi$ [rad]'); 
			h0.set_clim(vmin=-1.0,vmax=1.0); fig.colorbar(h0,label=r'$B_R$ [T]')
			h1=ax[0,1].pcolormesh(x,y,np.squeeze(beam_data.B_Z[:,:,j]).T,cmap='jet',shading='gouraud')
			ax[0,1].set_xlabel('R [m]'); ax[0,1].set_ylabel(r'$\phi$ [rad]'); 
			h1.set_clim(vmin=-1.0,vmax=1.0); fig.colorbar(h1,label=r'$B_Z$ [T]')
			h2=ax[1,0].pcolormesh(x,y,np.squeeze(beam_data.B_PHI[:,:,j]).T,cmap='jet',shading='gouraud')
			ax[1,0].set_xlabel('R [m]'); ax[1,0].set_ylabel(r'$\phi$ [rad]'); 
			h2.set_clim(vmin=-7.0,vmax=7.0); fig.colorbar(h2,label=r'$B_\phi$ [T]')
			h3=ax[1,1].pcolormesh(x,y,np.squeeze(b[:,:,j]).T,cmap='jet',shading='gouraud')
			ax[1,1].set_xlabel('R [m]'); ax[1,1].set_ylabel(r'$\phi$ [rad]'); 
			h3.set_clim(vmin=0.0,vmax=10.0); fig.colorbar(h3,label=r'$|B|$ [T]')
			if not args.lbackground:pyplot.show()
			if (args.lsave): fig.savefig(f'brphi_{args.brphi_index_phi:03d}_{args.beams3d_ext}.png', dpi=fig.dpi)
	sys.exit(0)