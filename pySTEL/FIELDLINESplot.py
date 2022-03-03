#!/usr/bin/env python3
import sys, os
os.environ['ETS_TOOLKIT'] = 'qt4'
import matplotlib
matplotlib.use("Qt4Agg")
import matplotlib.pyplot as _plt
import numpy as np                    #For Arrays
import numpy.matlib
from math import pi
#QT4
from PyQt4 import uic, QtGui
from PyQt4.QtGui import QMainWindow, QApplication, qApp, QVBoxLayout, QSizePolicy,QIcon
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
#QT5
#from PyQt5 import uic, QtGui, QtWidgets
#from PyQt5.QtWidgets import QMainWindow, QApplication, QVBoxLayout, QSizePolicy
#from PyQt5.QtGui import QIcon
#from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from libstell.fieldlines import read_fieldlines, calc_iota, calc_reff
from matplotlib.figure import Figure
from mpl_toolkits import mplot3d

try:
	qtCreatorPath=os.environ["STELLOPT_PATH"]
except KeyError:
	print("Please set environment variable STELLOPT_PATH")
	sys.exit(1)

qtCreatorFile = qtCreatorPath+"/pySTEL/FIELDLINESplot.ui" # Enter file here.
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

class MyApp(QMainWindow):
	def __init__(self):
		super(MyApp, self).__init__()
		self.ui = Ui_MainWindow()
		self.ui.setupUi(self) 
		self.setStyleSheet("background-color: white;");
		#self.ui.PlotButtons.setStyleSheet("background-color: white;");
		self.statusBar().showMessage('Ready')
		self.ui.plot_list = ['Summary', '-----1D-----', 'Iota', 'q',\
		'-----3D------', 'B_R', 'B_PHI', 'B_Z',\
		'---Special---', 'Poincaré', 'Poincaré |B|', 'Poincaré a',\
		'Poincaré Length']
		files = sorted(os.listdir('.'))
		for name in files:
			if(name[0:10]=='fieldlines'):
				self.ui.FileName.addItem(name)
		# Init
		self.fieldlines_data=read_fieldlines(self.ui.FileName.currentText())
		self.ui.PlotList.addItems(self.ui.plot_list)
		self.ui.PlotButtons.setEnabled(0)
		self.nlines = self.fieldlines_data['nlines']
		self.npoinc = self.fieldlines_data['npoinc']
		self.nsteps = self.fieldlines_data['nsteps']
		self.nr     = self.fieldlines_data['nr']
		self.nphi   = self.fieldlines_data['nphi']
		self.nz     = self.fieldlines_data['nz']
		self.ui.rhoslider.setMaximum(self.nr-1)
		self.ui.uslider.setMaximum(self.nphi-1)
		self.ui.vslider.setMaximum(self.nz-1)
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
		self.ui.savebutton.clicked.connect(self.plot_to_file)

	def FileSelect(self,i):
		self.fieldlines_data=read_fieldlines(self.ui.FileName.currentText())
		#self.ui.PlotList.addItems(self.ui.plot_list)
		self.nlines = self.fieldlines_data['nlines']
		self.npoinc = self.fieldlines_data['npoinc']
		self.nsteps = self.fieldlines_data['nsteps']
		self.nr     = self.fieldlines_data['nr']
		self.nphi   = self.fieldlines_data['nphi']
		self.nz     = self.fieldlines_data['nz']
		self.ui.rhoslider.setMaximum(self.nr-1)
		self.ui.uslider.setMaximum(self.nphi-1)
		self.ui.vslider.setMaximum(self.nz-1)
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
		elif (i < 4):
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
		elif (i>8):
			self.ui.poltor_button.setChecked(1)
			self.ui.rho_button.setEnabled(0)
			self.ui.pol_button.setEnabled(0)
			self.ui.tor_button.setEnabled(0)
			self.ui.flux_button.setEnabled(0)
			self.ui.poltor_button.setEnabled(0)
			self.ui.RZ_button.setEnabled(0)
			self.ui.ThreeD_button.setEnabled(0)
			self.ui.rhoslider.setEnabled(1)
			self.ui.uslider.setEnabled(1)
			self.ui.vslider.setEnabled(0)
			self.ui.rhoslider.setMaximum(self.nlines-1)
			self.ui.uslider.setMaximum(self.npoinc-1)
			#self.ui.vslider.setMaximum(self.nz-1)
			self.s=0; self.u=0; self.v=0;
		else:
			self.ui.rho_button.setChecked(1)
			self.CutSelect(self)
		self.update_plot(self)

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
		#self.ui.plot_widget.addWidget(self.canvas)
		plot_name = self.ui.PlotList.currentText();
		self.fig.clf()
		#self.fig.delaxes(self.ax)
		self.ax = self.fig.add_subplot(111)
		if (plot_name == 'Summary'):
			print(plot_name)
		elif (plot_name == 'Iota'):
			temp = calc_iota(self.fieldlines_data)
			#self.ax.errorbar(temp['rho'],temp['iota'],yerr=temp['iota_err'])
			self.ax.plot(temp['rho'],temp['iota'],'.k')
			self.ax.set_ylim(min(temp['iota']),max(temp['iota']))
			self.ax.set_xlim(0.0,1.25)
			self.ax.set_xlabel('Normalized Flux')
			self.ax.set_ylabel('iota')
			self.ax.set_title('Rotational Transform')
			#self.ax.set(xlabel='s',ylabel='iota',aspect='square')
		elif (plot_name == 'q'):
			temp = calc_iota(self.fieldlines_data)
			self.ax.plot(temp['rho'],1.0/temp['iota'],'.k')
			self.ax.set_xlabel('Normalized Flux')
			self.ax.set_ylabel('q')
			self.ax.set_title('Safety Factor')
		elif (plot_name == 'Poincaré'):
			k = self.u
			rmin = np.amin(self.fieldlines_data['raxis'])
			rmax = np.amax(self.fieldlines_data['raxis'])
			self.ax.plot(self.fieldlines_data['R_lines'][k:self.nsteps-1:self.npoinc,:],\
				self.fieldlines_data['Z_lines'][k:self.nsteps-1:self.npoinc,:],\
				'.k',markersize=0.1)
			self.ax.set_xlabel('R [m]')
			self.ax.set_ylabel('Z [m]')
			self.ax.set_title('Poincaré Plot')
			self.ax.set_aspect('equal')
			self.ax.set_xlim(rmin,rmax)
		elif (plot_name == 'Poincaré |B|'):
			k = self.u
			rmin = np.amin(self.fieldlines_data['raxis'])
			rmax = np.amax(self.fieldlines_data['raxis'])
			cmin = np.amin(self.fieldlines_data['B_lines'][k:self.npoinc,0])
			cmax = np.amax(self.fieldlines_data['B_lines'][k:self.npoinc,0])
			#cmin = cmin - 0.5*np.abs(cmin)
			#cmax = cmax + 0.5*np.abs(cmin)
			print(self.fieldlines_data['B_lines'][0,:])
			cax  = self.ax.scatter(self.fieldlines_data['R_lines'][k:self.nsteps-1:self.npoinc,:],\
				self.fieldlines_data['Z_lines'][k:self.nsteps-1:self.npoinc,:],\
				c=self.fieldlines_data['B_lines'][k:self.nsteps-1:self.npoinc,:], marker='.',s=0.3, \
				cmap='jet')
			self.ax.set_xlabel('R [m]')
			self.ax.set_ylabel('Z [m]')
			self.ax.set_title('Poincaré Plot')
			self.ax.set_aspect('equal')
			self.ax.set_xlim(rmin,rmax)
			self.fig.colorbar(cax)
			#cax.set_clim(cmin,cmax)
		elif (plot_name == 'Poincaré Length'):
			rho = self.fieldlines_data['L_lines']
			k = self.u
			rmin = np.amin(self.fieldlines_data['raxis'])
			rmax = np.amax(self.fieldlines_data['raxis'])
			cmin = np.amin(rho)
			cmax = np.amax(rho)
			n = np.rint((self.nsteps-1-k)/self.npoinc)
			rho2d = np.matlib.repmat(rho,int(n),1)
			cax  = self.ax.scatter(self.fieldlines_data['R_lines'][k:self.nsteps-1:self.npoinc,:],\
				self.fieldlines_data['Z_lines'][k:self.nsteps-1:self.npoinc,:],\
				c=rho2d, marker='.',s=0.3, \
				cmap='jet')
			self.ax.set_xlabel('R [m]')
			self.ax.set_ylabel('Z [m]')
			self.ax.set_title('Poincaré Plot ($L_{line}$)')
			self.ax.set_aspect('equal')
			self.ax.set_xlim(rmin,rmax)
			self.fig.colorbar(cax)
			cax.set_clim(cmin,cmax)
		elif (plot_name == 'Poincaré a'):
			rho = calc_reff(self.fieldlines_data)
			k = self.u
			rmin = np.amin(self.fieldlines_data['raxis'])
			rmax = np.amax(self.fieldlines_data['raxis'])
			cmin = np.amin(rho)
			cmax = np.amax(rho)
			n = np.rint((self.nsteps-1-k)/self.npoinc)
			rho2d = np.matlib.repmat(rho,int(n),1)
			cax  = self.ax.scatter(self.fieldlines_data['R_lines'][k:self.nsteps-1:self.npoinc,:],\
				self.fieldlines_data['Z_lines'][k:self.nsteps-1:self.npoinc,:],\
				c=rho2d, marker='.',s=0.3, \
				cmap='jet')
			self.ax.set_xlabel('R [m]')
			self.ax.set_ylabel('Z [m]')
			self.ax.set_title('Poincaré Plot ($a_{minor}$)')
			self.ax.set_aspect('equal')
			self.ax.set_xlim(rmin,rmax)
			self.fig.colorbar(cax)
			cax.set_clim(cmin,cmax)
		elif (plot_name[0] == '-'):
			print(plot_name)
		else:
			# First load the value based on plot
			if (plot_name=='B_R'):
				val = self.fieldlines_data['B_R']
			elif (plot_name=='B_PHI'):
				val = self.fieldlines_data['B_PHI']
			elif (plot_name=='B_Z'):
				val = self.fieldlines_data['B_Z']
			r = self.fieldlines_data['raxis']
			p = self.fieldlines_data['phiaxis']
			z = self.fieldlines_data['zaxis']
			nr2 = np.int(self.nr/2)
			nphi2 = np.int(self.nphi/2)
			nz2 = np.int(self.nz/2)
			cmin = self.fieldlines_data['B_PHI'][nr2,nphi2,nz2]
			cmax = cmin + np.abs(cmin)*0.75
			cmin = cmin - np.abs(cmin)*0.75
			# Now handle the type of plot
			if (self.ui.rho_button.isChecked()):
				self.ax.plot(r,val[:,self.u,self.v])
				self.ax.set_xlabel('R [m]')
			elif (self.ui.pol_button.isChecked()):
				self.ax.plot(p,val[self.s,:,self.v])
				self.ax.set_xlabel('Toridal Angle [rad]')
			elif (self.ui.tor_button.isChecked()):
				self.ax.plot(z,val[self.s,self.u,:])
				self.ax.set_xlabel('Z [m]')
			elif (self.ui.flux_button.isChecked()):
				cax=self.ax.pcolormesh(np.squeeze(val[self.s,:,:]),cmap='jet')
				self.fig.colorbar(cax)
				self.ax.set_ylabel('Toroidal Angle [rad]')
				self.ax.set_xlabel('Z [m]')
				cax.set_clim(cmin,cmax)
			elif (self.ui.poltor_button.isChecked()):
				cax=self.ax.pcolormesh(np.squeeze(val[:,self.u,:]),cmap='jet')
				self.fig.colorbar(cax)
				self.ax.set_xlabel('R [m]')
				self.ax.set_ylabel('Z [m]')
				cax.set_clim(cmin,cmax)
			elif (self.ui.RZ_button.isChecked()):
				cax = self.ax.pcolor(val[:,:,self.v],cmap='jet')
				self.fig.colorbar(cax)
				self.ax.set_xlabel('R [m]')
				self.ax.set_ylabel('Toroidal Angle [rad]')
				cax.set_clim(cmin,cmax)
			elif (self.ui.ThreeD_button.isChecked()):
				self.fig.delaxes(self.ax)
		self.canvas.draw()

	def plot_to_file(self,i):
		text = self.ui.saveas_filename.toPlainText();
		self.fig.savefig('./'+text, dpi=300)


if __name__ == "__main__":
	app = QApplication(sys.argv) 
	window = MyApp() 
	window.show() 
	sys.exit(app.exec_())
