#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Main routine
if __name__=="__main__":
	import sys
	from argparse import ArgumentParser
	import numpy as np
	import matplotlib.pyplot as pyplot
	from libstell.gist import GIST
	parser = ArgumentParser(description= 
		'''Provides a simple tool for reading and plotting GIST files.''')
	parser.add_argument("-g", "--gist", dest="gist_ext",
		help="GIST File name.", default = None)
	parser.add_argument("-p", "--plot", dest="lplot", action='store_true',
		help="Plot the GIST file.", default = False)
	parser.add_argument("-pp", "--plotproxy", dest="lplotproxy", action='store_true',
		help="Plot the Proxies.", default = False)
	gist_data = GIST()
	args = parser.parse_args()
	if args.gist_ext:
		try:
			gist_data.read_gist(args.gist_ext)
		except:
			print(f'Could not open GIST file: {args.gist_ext}')
			sys.exit(-1)
		if args.lplot:
			px = 1/pyplot.rcParams['figure.dpi']
			fig,axes = pyplot.subplots(2,4,figsize=(1024*px,768*px))
			pyplot.subplots_adjust(hspace=0.5,wspace=0.3)
			theta = np.linspace(-np.pi,np.pi,gist_data.gridpoints)
			# title
			fig.text(0.5, 0.925,rf'GIST (s={gist_data.s0},$\alpha$={gist_data.alpha0})',
				horizontalalignment='center',fontsize=24.0)
			# Figure 1
			axes[0,0].plot(theta,gist_data.g11,'k')
			axes[0,0].set_xlabel(r'Theta $\theta$')
			axes[0,0].set_title(r'$g_{11}$')
			axes[0,0].set_xlim([-np.pi,np.pi])
			# Figure 2
			axes[0,1].plot(theta,gist_data.g12,'k')
			axes[0,1].set_xlabel(r'Theta $\theta$')
			axes[0,1].set_title(r'$g_{12}$')
			axes[0,1].set_xlim([-np.pi,np.pi])
			# Figure 3
			axes[0,2].plot(theta,gist_data.g22,'k')
			axes[0,2].set_xlabel(r'Theta $\theta$')
			axes[0,2].set_title(r'$g_{22}$')
			axes[0,2].set_xlim([-np.pi,np.pi])
			# Figure 4
			axes[0,3].plot(theta,gist_data.Bhat,'k')
			axes[0,3].set_xlabel(r'Theta $\theta$')
			axes[0,3].set_title(r'$\hat{B}$')
			axes[0,3].set_xlim([-np.pi,np.pi])
			# Figure 5
			axes[1,0].plot(theta,gist_data.abs_jac,'k')
			axes[1,0].set_xlabel(r'Theta $\theta$')
			axes[1,0].set_title(r'$|\sqrt{g}|$')
			axes[1,0].set_xlim([-np.pi,np.pi])
			# Figure 6
			axes[1,1].plot(theta,gist_data.L1,'k')
			axes[1,1].set_xlabel(r'Theta $\theta$')
			axes[1,1].set_title(r'L1')
			axes[1,1].set_xlim([-np.pi,np.pi])
			# Figure 7
			axes[1,2].plot(theta,gist_data.L2,'k')
			axes[1,2].set_xlabel(r'Theta $\theta$')
			axes[1,2].set_title(r'L2')
			axes[1,2].set_xlim([-np.pi,np.pi])
			# Figure 8
			axes[1,3].plot(theta,gist_data.dBdt,'k')
			axes[1,3].set_xlabel(r'Theta $\theta$')
			axes[1,3].set_title(r'$dB/d\theta$')
			axes[1,3].set_xlim([-np.pi,np.pi])
			pyplot.show()
		if args.lplotproxy:
			px = 1/pyplot.rcParams['figure.dpi']
			fig,axes = pyplot.subplots(2,4,figsize=(1024*px,768*px))
			pyplot.subplots_adjust(hspace=0.3,wspace=0.4)
			theta = np.linspace(-np.pi,np.pi,gist_data.gridpoints)
			# title
			fig.text(0.5, 0.925,rf'Proxies (s={gist_data.s0},$\alpha$={gist_data.alpha0})',
				horizontalalignment='center',fontsize=24.0)
			# Figure 1
			axes[0,0].plot(theta,gist_data.calcProx1(),'k')
			axes[0,0].set_xlabel(r'Theta $\theta$')
			axes[0,0].set_title(r'Proxy1')
			axes[0,0].set_xlim([-np.pi,np.pi])
			# Figure 2
			axes[0,1].plot(theta,gist_data.calcProx1b(),'k')
			axes[0,1].set_xlabel(r'Theta $\theta$')
			axes[0,1].set_title(r'Proxy1b')
			axes[0,1].set_xlim([-np.pi,np.pi])
			# Figure 3
			axes[0,2].plot(theta,gist_data.calcProx1d(),'k')
			axes[0,2].set_xlabel(r'Theta $\theta$')
			axes[0,2].set_title(r'Proxy1d')
			axes[0,2].set_xlim([-np.pi,np.pi])
			# Figure 4
			axes[0,3].plot(theta,gist_data.calcProxL2(),'k')
			axes[0,3].set_xlabel(r'Theta $\theta$')
			axes[0,3].set_title(r'ProxyL2')
			axes[0,3].set_xlim([-np.pi,np.pi])
			# Remove axes
			gs = axes[1, 0].get_gridspec()
			axes[1,0].remove()
			axes[1,1].remove()
			axes[1,2].remove()
			axes[1,3].remove()
			axbig2 = fig.add_subplot(gs[1, 0:2])
			axbig1 = fig.add_subplot(gs[1, 2:4])
			# Figure 6
			axbig1.plot(theta,gist_data.Bhat,'b',label=r'$\hat{B}$')
			ax2=axbig1.twinx()
			ax2.plot(theta,gist_data.L2,'r',label='L2')
			axbig1.set_xlabel(r'Theta $\theta$')
			axbig1.set_ylabel(r'$\hat{B}$')
			axbig1.set_title(r'TEM Proxy Data')
			ax2.set_ylabel('L2')
			axbig1.set_xlim([-np.pi,np.pi])
			fig.text(0.02,0.05,rf'Bounce Proxy: {gist_data.calcProxTEMBounce():7.5f}')
			fig.text(0.02,0.02,rf'   Tau Proxy: {gist_data.calcProxTEMTau():7.5f}')
			# Figure 5
			axbig2.plot(theta,gist_data.calcProxTEMOverlap(),'k')
			axbig2.set_xlabel(r'Theta $\theta$')
			axbig2.set_title(r'Proxy TEM Overlap')
			axbig2.set_xlim([-np.pi,np.pi])
			pyplot.show()
	sys.exit(0)