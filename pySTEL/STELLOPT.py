#!/usr/bin/env python3
import sys, os
import matplotlib
matplotlib.use("Qt4Agg")
import matplotlib.pyplot as _plt
import numpy as np                    #For Arrays
from math import pi
from PyQt4 import uic, QtGui
from PyQt4.QtGui import QMainWindow, QApplication, qApp, QApplication, QVBoxLayout, \
                        QSizePolicy, QWidget, QFileDialog
from PyQt4.QtGui import QIcon, QTableWidget, QTableWidgetItem
from libstell.libstell import safe_open, read_indata_namelist, pmass, pcurr, piota, \
                              set_module_var, safe_close, cfunct, sfunct, isotoro, \
                              write_indata_namelist
from libstell.stellopt import read_stellopt_namelist
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits import mplot3d
# MayaVi stuff
#os.environ['ETS_TOOLKIT'] = 'qt4'
#from mayavi.core.ui.api import MayaviScene
#from mayavi import mlab

try:
	qtCreatorPath=os.environ["STELLOPT_PATH"]
except KeyError:
	print("Please set environment variable STELLOPT_PATH")
	sys.exit(1)

qtCreatorFile = qtCreatorPath+"/pySTEL/STELLOPT.ui" # Enter file here.
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

class MyApp(QMainWindow):
	def __init__(self):
		super(MyApp, self).__init__()
		self.ui = Ui_MainWindow()
		self.ui.setupUi(self) 
		self.setStyleSheet("background-color: white;")
		# Setup Defaults
		iunit = 55
		istat = 0
		self.indata=read_indata_namelist(iunit,istat) # dummy just to get going
		# Set the OPTIMUM DEFAULTS (will repalce with something like read_indata_namelist)
		self.optimum=read_stellopt_namelist(iunit,istat)
		# Setup Components
		self.ui.TableArrays.setRowCount(1)
		self.ui.TableArrays.setColumnCount(20)
		for num,item in enumerate(self.indata['am'], start=0):
			self.ui.TableArrays.setItem(0,num, QTableWidgetItem(str(item)))
		# Setup Plot
		self.fig = Figure(figsize=(2,2),dpi=100)
		self.ax = self.fig.add_subplot(111)
		self.canvas = FigureCanvas(self.fig)
		self.ui.Plotbox.addWidget(self.canvas)
		self.DrawArrays()
		# Callbacks (VMEC Tab)
		self.ui.TextMpol.editingFinished.connect(self.UpdateMpol)
		self.ui.TextNtor.editingFinished.connect(self.UpdateNtor)
		self.ui.ButtonLoadIndata.clicked.connect(self.LoadIndata)
		self.ui.ComboBoxArrays.currentIndexChanged.connect(self.UpdateArrays)
		self.ui.ComboBoxPType.currentIndexChanged.connect(self.UpdatePType)
		self.ui.TableArrays.cellChanged.connect(self.DataArrays)
		self.ui.TableNsArray.cellChanged.connect(self.NSArrays)
		self.ui.ButtonWriteIndata.clicked.connect(self.WriteIndata)
		# Callbacks (STELLOPT Tab)
		self.ui.ComboBoxOPTtype.currentIndexChanged.connect(self.UpdateOPTtype)

	def UpdateMpol(self):
		strtmp = self.ui.TextMpol.text()
		inttmp = int(strtmp)
		self.indata['mpol'] = inttmp
		set_module_var('vmec_input','mpol',inttmp)
		return

	def UpdateNtor(self):
		strtmp = self.ui.TextNtor.text()
		inttmp = int(strtmp)
		self.indata['ntor'] = inttmp
		set_module_var('vmec_input','ntor',inttmp)
		return

	def LoadIndata(self):
		# Handles loading an indata file.
		w = QWidget()
		w.resize(320, 240)
		w.setWindowTitle("Hello World!")
		filename = QFileDialog.getOpenFileName(w, 'Open File', '.')
		w.destroy
		# Now read the file
		iunit = 27
		istat = 0
		recl  = 1
		temp=safe_open(iunit,istat,filename,'old','formatted',recl,'sequential','none')
		self.indata=read_indata_namelist(iunit,istat)
		safe_close(iunit)
		# Now update the UI
		self.ui.TextTcon0.setText(str(self.indata['tcon0'])) 
		self.ui.TextDelt.setText(str(self.indata['delt'])) 
		self.ui.TextMpol.setText(str(self.indata['mpol'])) 
		self.ui.TextNtor.setText(str(self.indata['ntor'])) 
		self.ui.TextNtheta.setText(str(self.indata['ntheta'])) 
		self.ui.TextNzeta.setText(str(self.indata['nzeta'])) 
		self.ui.TextNfp.setText(str(self.indata['nfp'])) 
		self.ui.TextNvacSkip.setText(str(self.indata['nvacskip'])) 
		self.ui.TextMgridFile.setText(self.indata['mgrid_file']) 
		if self.indata['lfreeb']:
			self.ui.CheckboxLfreeb.setChecked(True)
		else:
			self.ui.CheckboxLfreeb.setChecked(False)
		self.ui.TextGamma.setText(str(self.indata['gamma']))
		self.ui.TextPhiEdge.setText(str(self.indata['phiedge']))
		self.ui.TextCurtor.setText(str(self.indata['curtor']))
		self.ui.TextPresScale.setText(str(self.indata['pres_scale'])) 
		self.ui.TextBloat.setText(str(self.indata['bloat'])) 
		self.ui.TextSpresPed.setText(str(self.indata['spres_ped']))
		self.ui.ComboBoxNcurr.setCurrentIndex(self.indata['ncurr'])
		# Handle the ComboBoxPType
		# Handle NS Array table widget
		ns_mask=self.indata['ns_array']!=0
		ns_array = self.indata['ns_array'][ns_mask]
		ftol_array = self.indata['ftol_array'][ns_mask]
		niter_array = self.indata['niter_array'][ns_mask]
		self.ui.TableNsArray.setColumnCount(len(ns_array))
		for num,item in enumerate(ns_array, start=0):
			self.ui.TableNsArray.setItem(0,num, QTableWidgetItem(str(item)))
			self.ui.TableNsArray.setItem(1,num, QTableWidgetItem(str(niter_array[num])))
			self.ui.TableNsArray.setItem(2,num, QTableWidgetItem(str(ftol_array[num])))
		self.ui.TableNsArray.show()
		self.UpdateArrays()

	def WriteIndata(self):
		# Handles loading an indata file.
		w = QWidget()
		w.resize(320, 240)
		w.setWindowTitle("Hello World!")
		filename = QFileDialog.getSaveFileName(w, 'Open File', '.')
		w.destroy
		# Update the module with all the not-updated values
		self.indata['tcon0'] = float(self.ui.TextTcon0.text())
		set_module_var('vmec_input','tcon0',self.indata['tcon0'])
		self.indata['delt'] = float(self.ui.TextDelt.text())
		set_module_var('vmec_input','delt',self.indata['delt'])
		self.indata['ntheta'] = int(self.ui.TextNtheta.text())
		set_module_var('vmec_input','ntheta',self.indata['ntheta'])
		self.indata['nzeta'] = int(self.ui.TextNzeta.text())
		set_module_var('vmec_input','nzeta',self.indata['nzeta'])
		self.indata['nfp'] = int(self.ui.TextNfp.text())
		set_module_var('vmec_input','nfp',self.indata['nfp'])
		self.indata['nvacskip'] = int(self.ui.TextNvacSkip.text())
		set_module_var('vmec_input','nvacskip',self.indata['nvacskip'])
		self.indata['gamma'] = float(self.ui.TextGamma.text())
		set_module_var('vmec_input','gamma',self.indata['gamma'])
		self.indata['phiedge'] = float(self.ui.TextPhiEdge.text())
		set_module_var('vmec_input','phiedge',self.indata['phiedge'])
		self.indata['pres_scale'] = float(self.ui.TextPresScale.text())
		set_module_var('vmec_input','pres_scale',self.indata['pres_scale'])
		self.indata['bloat'] = float(self.ui.TextBloat.text())
		set_module_var('vmec_input','bloat',self.indata['bloat'])
		self.indata['spres_ped'] = float(self.ui.TextSpresPed.text())
		set_module_var('vmec_input','spres_ped',self.indata['spres_ped'])
		self.indata['curtor'] = float(self.ui.TextCurtor.text())
		set_module_var('vmec_input','curtor',self.indata['curtor'])
		# Now write the file
		iunit = 27
		istat = 0
		recl  = 1000
		temp=safe_open(iunit,istat,filename,'unknown','formatted',recl,'sequential','none')
		self.indata=write_indata_namelist(iunit,istat) # Not working yet
		safe_close(iunit)

	def UpdatePType(self):
		data_name=self.ui.ComboBoxArrays.currentText()
		data_name = data_name.lower()
		type_name=self.ui.ComboBoxPType.currentText()
		type_name = type_name.lower()
		if data_name == 'am' or data_name == 'am_aux':
			set_module_var('vmec_input','pmass_type','                    ')
			set_module_var('vmec_input','pmass_type',type_name)
			self.indata['pmass_type'] = type_name
		elif data_name == 'ac' or data_name == 'ac_aux':
			set_module_var('vmec_input','pcurr_type','                    ')
			set_module_var('vmec_input','pcurr_type',type_name)
			self.indata['pcurr_type'] = type_name
		elif data_name == 'ai' or data_name == 'ai_aux':
			set_module_var('vmec_input','piota_type','                    ')
			set_module_var('vmec_input','piota_type',type_name)
			self.indata['piota_type'] = type_name
		else:
			return
		self.UpdateArrays()
		self.DrawArrays()

	def NSArrays(self):
		# Handles changes in NS matrix
		# Get changed item
		item = self.ui.TableNsArray.currentItem()
		if item is None:
			return
		col = self.ui.TableNsArray.currentColumn()
		row = self.ui.TableNsArray.currentRow()
		if row == 0: # NS
			self.indata['ns_array'][col]=int(item.text())
			set_module_var('vmec_input','ns_array',self.indata['ns_array'][:])
		elif row == 1: #NITER_ARRAY
			self.indata['niter_array'][col]=int(item.text())
			set_module_var('vmec_input','niter_array',self.indata['niter_array'][:])
		elif row ==2: #FTOL
			self.indata['ftol_array'][col]=float(item.text())
			set_module_var('vmec_input','ftol_array',self.indata['ftol_array'][:])

	def UpdateArrays(self):
		# Updates values in array box
		dex=self.ui.ComboBoxArrays.currentIndex()
		data_name=self.ui.ComboBoxArrays.currentText()
		data_name = data_name.lower()
		# Update PType accordingly
		strtmp = 'none'
		if data_name == 'am':
			strtmp = str(self.indata['pmass_type']).strip()
		if data_name == 'ai':
			strtmp = str(self.indata['piota_type']).strip()
		if data_name == 'ac' or data_name == 'ac_aux':
			strtmp = str(self.indata['pcurr_type']).strip()
		if data_name == 'am_aux':
			self.indata['pmass_type'] = 'akima_spline'
			set_module_var('vmec_input','pmass_type','akima_spline')
			strtmp = 'akima_spline'
		if data_name == 'ai_aux':
			self.indata['piota_type'] = 'akima_spline'
			set_module_var('vmec_input','piota_type','akima_spline')
			strtmp = 'akima_spline'
		if strtmp != 'none':
			self.ui.ComboBoxPType.setCurrentIndex(self.ui.ComboBoxPType.findText(strtmp))
		# Update the table
		self.ui.TableArrays.setRowCount(1)
		self.ui.TableArrays.setColumnCount(20)
		self.ui.TableArrays.setVerticalHeaderLabels('1')
		self.ui.TableArrays.setHorizontalHeaderLabels('0;1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20'.split(';'))
		if data_name in ['am','ac','ai']:
			for num,item in enumerate(self.indata[data_name], start=0):
				self.ui.TableArrays.setItem(0,num, QTableWidgetItem(str(item)))
		elif data_name in ['am_aux','ac_aux','ai_aux']:
			aux_mask = self.indata[data_name+'_s']>=0
			self.ui.TableArrays.setRowCount(2)
			self.ui.TableArrays.setColumnCount(50)
			if (np.count_nonzero(aux_mask) > 0):
				self.ui.TableArrays.setColumnCount(np.count_nonzero(aux_mask))
				for num,item in enumerate(self.indata[data_name+'_s'][aux_mask], start=0):
					self.ui.TableArrays.setItem(0,num, QTableWidgetItem(str(item)))
					self.ui.TableArrays.setItem(1,num, QTableWidgetItem(str(self.indata[data_name+'_f'][num])))
			else:
				self.ui.TableArrays.setColumnCount(6)
				self.ui.TableArrays.setItem(0,0, QTableWidgetItem('0.0'))
				self.ui.TableArrays.setItem(0,1, QTableWidgetItem('0.2'))
				self.ui.TableArrays.setItem(0,2, QTableWidgetItem('0.4'))
				self.ui.TableArrays.setItem(0,3, QTableWidgetItem('0.6'))
				self.ui.TableArrays.setItem(0,4, QTableWidgetItem('0.8'))
				self.ui.TableArrays.setItem(0,5, QTableWidgetItem('1.0'))
				self.ui.TableArrays.setItem(1,0, QTableWidgetItem('1.0'))
				self.ui.TableArrays.setItem(1,1, QTableWidgetItem('0.8'))
				self.ui.TableArrays.setItem(1,2, QTableWidgetItem('0.6'))
				self.ui.TableArrays.setItem(1,3, QTableWidgetItem('0.4'))
				self.ui.TableArrays.setItem(1,4, QTableWidgetItem('0.2'))
				self.ui.TableArrays.setItem(1,5, QTableWidgetItem('0.0'))
				self.indata[data_name+'_s']=np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
				self.indata[data_name+'_f']=np.array([1.0, 0.8, 0.6, 0.4, 0.2, 0.0])
				set_module_var('vmec_input',data_name+'_s',self.indata[data_name+'_s'])
				set_module_var('vmec_input',data_name+'_f',self.indata[data_name+'_f'])
		elif data_name == 'rbc' or data_name == 'zbs' or data_name == 'rbs' or data_name == 'zbc':

			self.ui.TableArrays.blockSignals(True)
			# Set Array Size
			self.ui.TableArrays.setRowCount(self.indata['ntor']*2+1)
			self.ui.TableArrays.setColumnCount(self.indata['mpol'])
			# Set Headings
			temp =''
			for i in range(2*self.indata['ntor']+1): temp = temp + str(i-self.indata['ntor'])+";"
			self.ui.TableArrays.setVerticalHeaderLabels(temp.split(";"))
			temp =''
			for i in range(self.indata['mpol']): temp = temp + str(i)+";"
			self.ui.TableArrays.setHorizontalHeaderLabels(temp.split(";"))
			# Set Values
			for i in range(2*self.indata['ntor']+1):
				for j in range(self.indata['mpol']):
					i2 = 101-self.indata['ntor'] + i
					val = str(self.indata[data_name][j,i2])
					self.ui.TableArrays.setItem(i,j, QTableWidgetItem(val))

			self.ui.TableArrays.blockSignals(False)
		self.ui.TableArrays.show()
		self.DrawArrays()

	def DataArrays(self):
		# Handles changes in matrix
		# Get type of array
		data_name = self.ui.ComboBoxArrays.currentText()
		data_name = data_name.lower()
		# Get changed item
		item = self.ui.TableArrays.currentItem()
		if item is None:
			return
		value = float(item.text())
		col = self.ui.TableArrays.currentColumn()
		row = self.ui.TableArrays.currentRow()
		if data_name in ['am_aux','ac_aux','ai_aux']:
			if row == 0:
				self.indata[data_name+'_s'][col]=value
				set_module_var('vmec_input',data_name+'_s',self.indata[data_name+'_s'][:])
			else:
				self.indata[data_name+'_f'][col]=value
				set_module_var('vmec_input',data_name+'_f',self.indata[data_name+'_f'][:])
		elif data_name in ['am','ac','ai']:
			self.indata[data_name][col]=value
			set_module_var('vmec_input',data_name,self.indata[data_name][:])
		elif data_name == 'rbc' or data_name == 'zbs' or data_name == 'rbs' or data_name == 'zbc':
			row2 = 101 - self.indata['ntor'] + row
			self.indata[data_name][col][row2]=value
		self.DrawArrays()

	def DrawArrays(self):
		# Determine type of array
		data_name = self.ui.ComboBoxArrays.currentText()
		data_name = data_name.lower()
		# Handle plots / values
		self.fig.clf()
		self.ax = self.fig.add_subplot(111)
		s = np.ndarray((99,1))
		f = np.ndarray((99,1))
		for i in range(99): s[i]=(i)/98.0
		if data_name[0:2] == 'am':
			for i,xx in enumerate(s):
				f[i] = pmass(xx)/(pi*4E-7)
				self.ax.set_ylabel('Pressure')
		elif data_name[0:2] == 'ac':
			for i,xx in enumerate(s):
				f[i] = pcurr(xx)
				self.ax.set_ylabel('Current')
		elif data_name[0:2] == 'ai':
			for i,xx in enumerate(s):
				f[i] = piota(xx)
				self.ax.set_ylabel('Iota')
		elif data_name == 'rbc' or data_name == 'zbs' or data_name == 'rbs' or data_name == 'zbc':
			mnmax = (2*self.indata['ntor']+1)*(self.indata['mpol'])
			nu = 4*self.indata['mpol']
			nv = 4*self.indata['ntor']
			if nu < 64: nu=64
			if nv < 32: nv=32
			mn = 0
			xn = np.ndarray((mnmax,1))
			xm = np.ndarray((mnmax,1))
			rmnc = np.ndarray((1,mnmax))
			zmns = np.ndarray((1,mnmax))
			rmns = np.ndarray((1,mnmax))
			zmnc = np.ndarray((1,mnmax))
			for i in range(2*self.indata['ntor']+1):
				for j in range(self.indata['mpol']):
					i2 = 101-self.indata['ntor'] + i
					xn[mn] = (i-self.indata['ntor'])*self.indata['nfp']
					xm[mn] = j
					rmnc[0,mn] = self.indata['rbc'][j,i2]
					zmns[0,mn] = self.indata['zbs'][j,i2]
					rmns[0,mn] = self.indata['rbs'][j,i2]
					zmnc[0,mn] = self.indata['zbc'][j,i2]
					mn = mn + 1
			theta = np.ndarray((nu,1))
			zeta = np.ndarray((nv,1))
			for j in range(nu): theta[j]=2*pi*j/(nu-1)
			for j in range(nv): zeta[j]=pi*j/(nv)/self.indata['nfp']
			r=cfunct(theta,zeta,rmnc,xm,xn)
			z=sfunct(theta,zeta,zmns,xm,xn)
			self.ax = isotoro(r,z,zeta,0,fig=self.fig)
			self.ax.grid(False)
			self.ax.set_axis_off()
			self.canvas.draw()
			return
		fmax = max(f)
		if fmax == 0:
			fmax = 1
		self.ax.plot(s,f)
		self.ax.set_xlabel('Norm. Tor. Flux (s)')
		self.ax.set_aspect('auto')
		if data_name[3:6] == 'aux':
			self.ax.plot(self.indata[data_name+'_s'],self.indata[data_name+'_f'],'o')
		self.canvas.draw()

	def UpdateOPTtype(self):
		# Determine type of Optimization
		self.optimum['OPT_TYPE'] = self.ui.ComboBoxOPTtype.currentText()
		self.optimum['OPT_TYPE'] = self.optimum['OPT_TYPE'].upper()
		if self.optimum['OPT_TYPE'] == 'ONE_ITER':
			fields = ['NOPTIMIZERS']
			self.ui.TableOPTtype.setRowCount(1)
		elif self.optimum['OPT_TYPE'] in ['LMDIF','LMDIF_BOUNDED']:
			self.ui.TableOPTtype.setRowCount(8)
			fields = ['NFUNC_MAX','FTOL','XTOL','GTOL','EPSFCN','FACTOR','MODE','NOPTIMIZERS']
		elif self.optimum['OPT_TYPE'] == 'GADE':
			self.ui.TableOPTtype.setRowCount(6)
			fields = ['NFUNC_MAX','FACTOR','CR_STRATEGY','MODE','NPOPULATION','NOPTIMIZERS']
		for i,item in enumerate(fields):
			self.ui.TableOPTtype.setItem(i,0, QTableWidgetItem(item))
			val = str(self.optimum[item])
			self.ui.TableOPTtype.setItem(i,1, QTableWidgetItem(val))










if __name__ == "__main__":
	app = QApplication(sys.argv) 
	window = MyApp() 
	window.show() 
	sys.exit(app.exec_())
