#!/usr/bin/env python3
import sys, os
import matplotlib
matplotlib.use("Qt4Agg")
import matplotlib.pyplot as _plt
import numpy as np                    #For Arrays
from math import pi
from PyQt4 import uic, QtGui
from PyQt4.QtGui import QMainWindow, QApplication, qApp, QApplication, QVBoxLayout, QSizePolicy
from PyQt4.QtGui import QIcon
from libstell.libstell import read_vmec, cfunct, sfunct, torocont, isotoro, calc_jll
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits import mplot3d

try:
	qtCreatorPath=os.environ["STELLOPT_PATH"]
except KeyError:
	print("Please set environment variable STELLOPT_PATH")
	sys.exit(1)

#qtCreatorFile = "/u/slazerso/src/STELLOPT_GCC/pySTEL/VMECplot.ui" # Enter file here.
qtCreatorFile = qtCreatorPath+"/pySTEL/VMECplot.ui" # Enter file here.
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

class MyApp(QMainWindow):
	def __init__(self):
		super(MyApp, self).__init__()
		self.ui = Ui_MainWindow()
		self.ui.setupUi(self) 
		#self.setStyleSheet("background-color: white;");
		self.statusBar().showMessage('Ready')
		self.ui.plot_list = ['Summary','-----1D-----','Iota','q','Pressure',\
		'<Buco>','<Bvco>','<jcuru>','<jcurv>','<j.B>',  '-----3D------','|B|','sqrt(g)',\
		'B^u','B^v','B_s','B_u','B_v','j^u','j^v', 'jll', 'j.B','---Special---','LPK']
		files = os.listdir('.')
		for name in files:
			if(name[0:4]=='wout'):
				self.ui.FileName.addItem(name)
		# Init
		self.vmec_data=read_vmec(self.ui.FileName.currentText())
		self.ui.PlotList.addItems(self.ui.plot_list)
		self.ui.PlotButtons.setEnabled(0)
		self.ns = self.vmec_data['ns']
		self.nu = self.vmec_data['mpol']*4
		self.nv = self.vmec_data['ntor']*4*self.vmec_data['nfp']
		self.nv2 = self.vmec_data['ntor']*4
		self.TransformVMEC(self)
		self.s=0
		self.u=0
		self.v=0
		self.ui.rhoslider.setMaximum(self.ns-1)
		self.ui.uslider.setMaximum(self.nu-1)
		self.ui.vslider.setMaximum((self.nv/self.vmec_data['nfp']))
		# Plot figure
		self.fig = Figure(figsize=(2,2),dpi=100)
		self.ax = self.fig.add_subplot(111)
		self.canvas = FigureCanvas(self.fig)
		self.ui.plot_widget.addWidget(self.canvas)
		#self.canvas.draw()
		# Callbacks		
		self.ui.FileName.currentIndexChanged.connect(self.FileSelect)
		self.ui.PlotList.currentIndexChanged.connect(self.PlotSelect)
		self.ui.rho_button.toggled.connect(self.CutSelect)
		self.ui.pol_button.toggled.connect(self.CutSelect)
		self.ui.tor_button.toggled.connect(self.CutSelect)
		self.ui.flux_button.toggled.connect(self.CutSelect)
		self.ui.poltor_button.toggled.connect(self.CutSelect)
		self.ui.RZ_button.toggled.connect(self.CutSelect)
		self.ui.ThreeD_button.toggled.connect(self.CutSelect)
		self.ui.rhoslider.valueChanged.connect(self.CutSelect)
		self.ui.uslider.valueChanged.connect(self.CutSelect)
		self.ui.vslider.valueChanged.connect(self.CutSelect)

	def FileSelect(self,i):
		self.vmec_data=read_vmec(self.ui.FileName.currentText())
		#self.ui.PlotList.addItems(self.ui.plot_list)
		self.ns = self.vmec_data['ns']
		self.nu = self.vmec_data['mpol']*4
		self.nv = self.vmec_data['ntor']*4
		self.TransformVMEC(self)
		self.s=0
		self.u=0
		self.v=0
		self.ui.rhoslider.setMaximum(self.ns-1)
		self.ui.uslider.setMaximum(self.nu-1)
		self.ui.vslider.setMaximum(self.nv-1)
		self.ui.PlotButtons.setEnabled(0)
		self.ui.PlotList.setCurrentIndex(0)
		self.update_plot(self)

	def PlotSelect(self,i):
		plot_name = self.ui.PlotList.currentText()
		# Handle Enable/Disable
		self.ui.PlotButtons.setEnabled(1)
		self.ui.rho_button.setEnabled(1)
		self.ui.pol_button.setEnabled(1)
		self.ui.tor_button.setEnabled(1)
		self.ui.flux_button.setEnabled(1)
		self.ui.poltor_button.setEnabled(1)
		self.ui.RZ_button.setEnabled(1)
		self.ui.ThreeD_button.setEnabled(1)
		self.ui.rhoslider.setEnabled(1)
		self.ui.uslider.setEnabled(1)
		self.ui.vslider.setEnabled(1)
		if (i==0):
			self.ui.PlotButtons.setEnabled(0)
			self.ui.rhoslider.setEnabled(0)
			self.ui.uslider.setEnabled(0)
			self.ui.vslider.setEnabled(0)
		elif (i<10):
			self.ui.rho_button.setChecked(1)
			self.ui.pol_button.setEnabled(0)
			self.ui.tor_button.setEnabled(0)
			self.ui.flux_button.setEnabled(0)
			self.ui.poltor_button.setEnabled(0)
			self.ui.RZ_button.setEnabled(0)
			self.ui.ThreeD_button.setEnabled(0)
			self.ui.rhoslider.setEnabled(0)
			self.ui.uslider.setEnabled(0)
			self.ui.vslider.setEnabled(0)
			self.s=0; self.u=0; self.v=0;
			self.update_plot(self)
		else:
			self.ui.rho_button.setChecked(1)
			self.CutSelect(self)

	def CutSelect(self,i):
		self.ui.rhoslider.setEnabled(1)
		self.ui.uslider.setEnabled(1)
		self.ui.vslider.setEnabled(1)
		#print(self.ui.rho_button.isChecked())
		if (self.ui.rho_button.isChecked()):
			self.ui.rhoslider.setEnabled(0)
			self.u = self.ui.uslider.value()
			self.v = self.ui.vslider.value()
		elif (self.ui.pol_button.isChecked()):
			self.ui.uslider.setEnabled(0)
			self.s = self.ui.rhoslider.value()
			self.v = self.ui.vslider.value()
		elif (self.ui.tor_button.isChecked()):
			self.ui.vslider.setEnabled(0)
			self.s = self.ui.rhoslider.value()
			self.u = self.ui.uslider.value()
		elif (self.ui.flux_button.isChecked()):
			self.ui.uslider.setEnabled(0)
			self.ui.vslider.setEnabled(0)
			self.s = self.ui.rhoslider.value()
		elif (self.ui.poltor_button.isChecked()):
			self.ui.rhoslider.setEnabled(0)
			self.ui.vslider.setEnabled(0)
			self.u = self.ui.uslider.value()
		elif (self.ui.RZ_button.isChecked()):
			self.ui.rhoslider.setEnabled(0)
			self.ui.uslider.setEnabled(0)
			self.v = self.ui.vslider.value()
		elif (self.ui.ThreeD_button.isChecked()):
			self.ui.uslider.setEnabled(0)
			self.ui.vslider.setEnabled(0)
			self.s = self.ui.rhoslider.value()
		self.update_plot(self)

	def update_plot(self,i):

		plot_name = self.ui.PlotList.currentText();
		self.fig.clf()
		#self.fig.delaxes(self.ax)
		self.ax = self.fig.add_subplot(111)
		if (plot_name == 'Summary'):
			print(plot_name)
		elif (plot_name == 'Iota'):
			self.ax.plot(self.nflux,self.vmec_data['iotaf'])
			self.ax.set_xlabel('Normalized Flux')
			self.ax.set_ylabel('iota')
			self.ax.set_title('Rotational Transform')
			#self.ax.set(xlabel='s',ylabel='iota',aspect='square')
		elif (plot_name == 'q'):
			self.ax.plot(self.nflux,1.0/self.vmec_data['iotaf'])
			self.ax.set_xlabel('Normalized Flux')
			self.ax.set_ylabel('q')
			self.ax.set_title('Safety Factor')
		elif (plot_name == 'Pressure'):
			self.ax.plot(self.nflux,self.vmec_data['presf']/1000)
			self.ax.set_xlabel('Normalized Flux')
			self.ax.set_ylabel('Pressure [kPa]')
			self.ax.set_title('Pressure Profile')
		elif (plot_name == '<Buco>'):
			self.ax.plot(self.nflux,self.vmec_data['buco'])
			self.ax.set_xlabel('Normalized Flux')
			self.ax.set_ylabel('<B^u> [T]')
			self.ax.set_title('Flux surface Averaged B^u')
		elif (plot_name == '<Bvco>'):
			self.ax.plot(self.nflux,self.vmec_data['bvco'])
			self.ax.set_xlabel('Normalized Flux')
			self.ax.set_ylabel('<B^v> [T]')
			self.ax.set_title('Flux surface Averaged B^v')
		elif (plot_name == '<jcuru>'):
			self.ax.plot(self.nflux,self.vmec_data['jcuru']/1000)
			self.ax.set_xlabel('Normalized Flux')
			self.ax.set_ylabel('<j^u> [kA/m^2]')
			self.ax.set_title('Flux surface Averaged j^u')
		elif (plot_name == '<jcurv>'):
			self.ax.plot(self.nflux,self.vmec_data['jcurv']/1000)
			self.ax.set_xlabel('Normalized Flux')
			self.ax.set_ylabel('<j^v> [kA/m^2]')
			self.ax.set_title('Flux surface Averaged j^v')
		elif (plot_name == '<j.B>'):
			self.ax.plot(self.nflux,self.vmec_data['jdotb']/1000)
			self.ax.set_xlabel('Normalized Flux')
			self.ax.set_ylabel('<j.B> [T*kA/m^2]')
			self.ax.set_title('Flux surface Averaged j.B')
		elif (plot_name == 'LPK'):
			self.ax.plot(self.r[self.ns-1,:,0],self.z[self.ns-1,:,0],color='red')
			self.ax.plot(self.r[0,0,0],self.z[0,0,0],'+',color='red')
			self.ax.plot(self.r[self.ns-1,:,int(self.nv2/4)],self.z[self.ns-1,:,int(self.nv2/4)],color='green')
			self.ax.plot(self.r[0,0,int(self.nv2/4)],self.z[0,0,int(self.nv2/4)],'+',color='green')
			self.ax.plot(self.r[self.ns-1,:,int(self.nv2/2)],self.z[self.ns-1,:,int(self.nv2/2)],color='blue')
			self.ax.plot(self.r[0,0,int(self.nv2/2)],self.z[0,0,int(self.nv2/2)],'+',color='blue')
			self.ax.set_xlabel('R [m]')
			self.ax.set_ylabel('Z [m]')
			self.ax.set_title('LPK Plot')
			self.ax.set_aspect('equal')
		elif (plot_name[0] == '-'):
			print(plot_name)
		else:
			# First load the value based on plot
			if (plot_name=='|B|'):
				val = self.b
			elif (plot_name=='sqrt(g)'):
				val = self.g
			elif (plot_name=='B^u'):
				val = self.bu
			elif (plot_name=='B^v'):
				val = self.bv
			elif (plot_name=='B_s'):
				val = self.b_s
			elif (plot_name=='B_u'):
				val = self.b_u
			elif (plot_name=='B_v'):
				val = self.b_v
			elif (plot_name=='j^u'):
				val = self.cu/self.g
			elif (plot_name=='j^v'):
				val = self.cv/self.g
			elif (plot_name=='jll'):
				val = calc_jll(self.vmec_data, self.theta, self.zeta)
			elif (plot_name=='j.B'):
				val = (self.cu*self.bu+self.cv*self.bv)/self.g
			# Now handle the type of plot
			if (self.ui.rho_button.isChecked()):
				self.ax.plot(self.nflux,val[:,self.u,self.v])
				self.ax.set_xlabel('Normalized Flux')
			elif (self.ui.pol_button.isChecked()):
				self.ax.plot(self.theta,val[self.s,:,self.v])
				self.ax.set_xlabel('Poloidal Angle [rad]')
			elif (self.ui.tor_button.isChecked()):
				self.ax.plot(self.zeta2,val[self.s,self.u,0:self.nv2+1])
				self.ax.set_xlabel('Toroidal Angle [rad]')
			elif (self.ui.flux_button.isChecked()):
				self.ax.pcolormesh(np.squeeze(val[self.s,:,0:self.nv2+1]),cmap='jet',shading='gouraud')
				self.ax.set_xlabel('Toroidal Angle [rad]')
				self.ax.set_ylabel('Poloidal Angle [rad]')
			elif (self.ui.poltor_button.isChecked()):
				self.ax.pcolormesh(np.squeeze(val[:,self.u,0:self.nv2+1]),cmap='jet',shading='gouraud')
				self.ax.set_xlabel('Toroidal Angle [rad]')
				self.ax.set_ylabel('Normalized Flux')
			elif (self.ui.RZ_button.isChecked()):
				#self.ax.pcolormesh(self.r[:,:,self.v],self.z[:,:,self.v],val[:,:,self.v],cmap='jet',shading='gouraud')
				cax = self.ax.pcolor(self.r[:,:,self.v],self.z[:,:,self.v],val[:,:,self.v],cmap='jet')
				self.fig.colorbar(cax)
				self.ax.set_xlabel('R [m]')
				self.ax.set_ylabel('Z [m]')
				self.ax.set_aspect('equal')
			elif (self.ui.ThreeD_button.isChecked()):
				self.fig.delaxes(self.ax)
				self.canvas.draw()
				self.ax = isotoro(self.r,self.z,self.zeta,self.s,val,fig=self.fig)
				self.ax.grid(False)
				self.ax.set_axis_off()
		self.canvas.draw()

	def TransformVMEC(self, i):
		self.nflux = np.ndarray((self.ns,1))
		self.theta = np.ndarray((self.nu,1))
		self.zeta = np.ndarray((self.nv,1))
		self.zeta2 = np.ndarray((self.nv2+1,1))
		for j in range(self.ns): self.nflux[j]=j/(self.ns-1)
		for j in range(self.nu): self.theta[j]=2*pi*j/(self.nu-1)
		for j in range(self.nv): self.zeta[j]=2*pi*j/((self.nv-1))
		self.zeta2=self.zeta[0:self.nv2+1]
		self.r=cfunct(self.theta,self.zeta,self.vmec_data['rmnc'],self.vmec_data['xm'],self.vmec_data['xn'])
		self.z=sfunct(self.theta,self.zeta,self.vmec_data['zmns'],self.vmec_data['xm'],self.vmec_data['xn'])
		self.b=cfunct(self.theta,self.zeta,self.vmec_data['bmnc'],self.vmec_data['xm'],self.vmec_data['xn'])
		self.g=cfunct(self.theta,self.zeta,self.vmec_data['gmnc'],self.vmec_data['xm'],self.vmec_data['xn'])
		self.bu=cfunct(self.theta,self.zeta,self.vmec_data['bsupumnc'],self.vmec_data['xm'],self.vmec_data['xn'])
		self.bv=cfunct(self.theta,self.zeta,self.vmec_data['bsupvmnc'],self.vmec_data['xm'],self.vmec_data['xn'])
		self.cu=cfunct(self.theta,self.zeta,self.vmec_data['currumnc'],self.vmec_data['xm'],self.vmec_data['xn'])
		self.cv=cfunct(self.theta,self.zeta,self.vmec_data['currvmnc'],self.vmec_data['xm'],self.vmec_data['xn'])
		self.b_s=sfunct(self.theta,self.zeta,self.vmec_data['bsubsmns'],self.vmec_data['xm'],self.vmec_data['xn'])
		self.b_u=cfunct(self.theta,self.zeta,self.vmec_data['bsubumnc'],self.vmec_data['xm'],self.vmec_data['xn'])
		self.b_v=cfunct(self.theta,self.zeta,self.vmec_data['bsubvmnc'],self.vmec_data['xm'],self.vmec_data['xn'])


if __name__ == "__main__":
	app = QApplication(sys.argv) 
	window = MyApp() 
	window.show() 
	sys.exit(app.exec_())
