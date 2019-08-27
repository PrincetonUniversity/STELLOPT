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
                              write_indata_namelist, read_vmec, cfunct, sfunct
from libstell.stellopt import read_stellopt_namelist, write_stellopt_namelist, read_stellopt
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
		self.indata['ntor']=4
		self.indata['mpol']=9
		self.indata['nfp']=3
		self.ui.tabMain.setTabEnabled(0,False)
		self.ui.tabMain.setTabEnabled(2,False)
		# Set the OPTIMUM DEFAULTS (will repalce with something like read_indata_namelist)
		self.optimum=read_stellopt_namelist(iunit,istat)
		# Setup Components
		self.ui.TableArrays.setRowCount(1)
		self.ui.TableArrays.setColumnCount(20)
		for num,item in enumerate(self.indata['am'], start=0):
			self.ui.TableArrays.setItem(0,num, QTableWidgetItem(str(item)))
		# Setup Plot VMEC
		self.fig = Figure(figsize=(2,2),dpi=100)
		self.ax = self.fig.add_subplot(111)
		self.canvas = FigureCanvas(self.fig)
		self.ui.Plotbox.addWidget(self.canvas)
		# Setup Plot STELLOPT
		self.fig2 = Figure(figsize=(2,2),dpi=100)
		self.ax2 = self.fig2.add_subplot(111)
		self.canvas2 = FigureCanvas(self.fig2)
		self.ui.OPTplot_box.addWidget(self.canvas2)
		# Setup STELLOPT Pannels
		self.UpdateOPTtype()
		self.UpdateOPTVarsScalar()
		self.UpdateOPTVarsProf()
		self.UpdateOPTVarsExtcur()
		# Callbacks (VMEC Tab)
		self.ui.TextMpol.editingFinished.connect(self.UpdateMpol)
		self.ui.TextNtor.editingFinished.connect(self.UpdateNtor)
		self.ui.TextNfp.editingFinished.connect(self.UpdateNfp)
		self.ui.ButtonLoadIndata.clicked.connect(self.LoadIndata)
		self.ui.ComboBoxArrays.currentIndexChanged.connect(self.UpdateArrays)
		self.ui.ComboBoxPType.currentIndexChanged.connect(self.UpdatePType)
		self.ui.TableArrays.cellChanged.connect(self.DataArrays)
		self.ui.TableNsArray.cellChanged.connect(self.NSArrays)
		self.ui.ButtonWriteIndata.clicked.connect(self.WriteIndata)
		# Callbacks (STELLOPT Tab)
		self.ui.ComboBoxOPTtype.currentIndexChanged.connect(self.UpdateOPTtype)
		self.ui.TableOPTtype.cellChanged.connect(self.OPTArrays)
		self.ui.tabStelVars.currentChanged.connect(self.StelVarsTab)
		self.ui.TableOPTVarsScalar.cellChanged.connect(self.OPTVarsScalar)
		self.ui.TableOPTVarsProf.cellChanged.connect(self.OPTVarsProf)
		self.ui.comboBoxStelVarsProfType.currentIndexChanged.connect(self.UpdateOPTVarsProf)
		self.ui.TableOPTVarsExtcur.cellChanged.connect(self.OPTVarsExtcur)
		self.ui.comboBoxStelVarsBoundType.currentIndexChanged.connect(self.UpdateOPTVarsBound)
		self.ui.TableOPTVarsBound.cellChanged.connect(self.OPTVarsBound)
		# Callbacks (OPT_plot Tab)
		self.ui.ButtonLoadSTELLOPT.clicked.connect(self.LoadSTELLOPT)
		self.ui.ComboBoxOPTplot_type.currentIndexChanged.connect(self.UpdateOptplot)

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
		return

	def UpdateNfp(self):
		strtmp = self.ui.TextNfp.text()
		inttmp = int(strtmp)
		self.indata['nfp'] = inttmp
		set_module_var('vmec_input','nfp',inttmp)
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
		self.indata['mgrid_file'] = self.ui.TextMgridFile.text()
		set_module_var('vmec_input','mgrid_file',self.indata['mgrid_file'])
		# Update NCURR
		dex = self.ui.ComboBoxNcurr.CurrentIndex()
		self.indata['ncurr'] = int(dex)
		set_module_var('vmec_input','ncurr',self.indata['ncurr'])
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
		self.fig.delaxes(self.ax)
		self.fig.clf()
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
			for j in range(nv): zeta[j]=pi*j/(nv*self.indata['nfp'])
			r=cfunct(theta,zeta,rmnc,xm,xn)
			z=sfunct(theta,zeta,zmns,xm,xn)
			self.ax = self.fig.add_subplot(111,projection='3d')
			self.ax = isotoro(r,z,zeta,0,fig=self.fig)
			self.ax.grid(False)
			self.ax.set_axis_off()
			self.canvas.draw()
			return
		fmax = max(f)
		if fmax == 0:
			fmax = 1
		self.ax = self.fig.add_axes([0, 0, 1, 1])
		#self.ax = self.fig.add_subplot(111)
		self.ax.plot(s,f)
		self.ax.set_xlabel('Norm. Tor. Flux (s)')
		self.ax.set_aspect('auto')
		if data_name[3:6] == 'aux':
			self.ax.plot(self.indata[data_name+'_s'],self.indata[data_name+'_f'],'o')
		self.canvas.draw()

	def UpdateOPTtype(self):
		self.ui.TableOPTtype.blockSignals(True)
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
		elif self.optimum['OPT_TYPE'] == 'PSO':
			self.ui.TableOPTtype.setRowCount(6)
			fields = ['NFUNC_MAX','FACTOR','CR_STRATEGY','MODE','NPOPULATION','NOPTIMIZERS']
		self.ui.TableOPTtype.setVerticalHeaderLabels(fields)
		for i,item in enumerate(fields):
			val = str(self.optimum[item])
			self.ui.TableOPTtype.setItem(i,0, QTableWidgetItem(val))
		self.ui.TableOPTtype.blockSignals(False)

	def OPTArrays(self):
		# Handles changes in OPT table
		# Get changed item
		item = self.ui.TableOPTtype.currentItem()
		if item is None:
			return
		col = self.ui.TableOPTtype.currentColumn()
		row = self.ui.TableOPTtype.currentRow()
		field = self.ui.TableOPTtype.verticalHeaderItem(row).text()
		print(field)
		if field == 'OPT_TYPE':
			self.optimum[field]=item.text()
		elif field in ['NOPTIMIZERS','NFUNC_MAX','MODE','CR_STRATEGY','NPOPULATION']:
			self.optimum[field]=int(item.text())
		else:
			self.optimum[field]=float(item.text())

	def StelVarsTab(self):
		# Handles updates to the tabs
		# Get page index
		item = self.ui.tabStelVars.currentIndex()
		print(item)
		if item == 0:
			self.UpdateOPTVarsScalar()
		elif item == 1:
			self.UpdateOPTVarsProf()
		elif item == 2:
			self.UpdateOPTVarsExtcur()
		elif item == 3:
			self.UpdateOPTVarsBound()


	def UpdateOPTVarsScalar(self):
		self.ui.TableOPTVarsScalar.blockSignals(True)
		#Update the Scalar Table
		for i,item in enumerate(['PHIEDGE','PRES_SCALE','CURTOR']):
			if self.optimum['L'+item+'_OPT'] == 1:
				self.ui.TableOPTVarsScalar.setItem(i,0,QTableWidgetItem('T'))
			else:
				self.ui.TableOPTVarsScalar.setItem(i,0,QTableWidgetItem('F'))
			self.ui.TableOPTVarsScalar.setItem(i,1,QTableWidgetItem(str(self.indata[item.lower()])))
			self.ui.TableOPTVarsScalar.setItem(i,2,QTableWidgetItem(str(self.optimum['D'+item+'_OPT'])))
			self.ui.TableOPTVarsScalar.setItem(i,3,QTableWidgetItem(str(self.optimum[item+'_MIN'])))
			self.ui.TableOPTVarsScalar.setItem(i,4,QTableWidgetItem(str(self.optimum[item+'_MAX'])))
		self.ui.TableOPTVarsScalar.blockSignals(False)

	def OPTVarsScalar(self):
		# Handles changes in Scalar table
		# Get changed item
		item = self.ui.TableOPTVarsScalar.currentItem()
		if item is None:
			return
		col = self.ui.TableOPTVarsScalar.currentColumn()
		row = self.ui.TableOPTVarsScalar.currentRow()
		field = self.ui.TableOPTVarsScalar.verticalHeaderItem(row).text()
		# Column determins value
		field = field.upper()
		if col == 0:
			field = 'L'+field+'_OPT'
			if item.text() == 'T':
				val = 1
			else:
				val = 0
		if col == 1:
			return
		elif col == 2:
			field = 'D'+field+'_OPT'
			val = float(item.text())
		elif col == 3:
			field = field+'_MIN'
			val = float(item.text())
		elif col == 4:
			field = field+'_MAX'
			val = float(item.text())
		print(field, val)
		self.optimum[field]=val

	def UpdateOPTVarsProf(self):
		self.ui.TableOPTVarsProf.blockSignals(True)
		# Figure out which array to deal with
		data_name = self.ui.comboBoxStelVarsProfType.currentText()
		if data_name == 'Pressure (AM)':
			field = 'AM'
		elif data_name == 'Current (AC)':
			field = 'AC'
		elif data_name == 'Iota (AI)':
			field = 'AI'
		if field in ['AM','AC','AI']:
			nrows = 11;
		self.ui.TableOPTVarsProf.setRowCount(nrows)
		for i in range(nrows):
			self.ui.TableOPTVarsProf.setVerticalHeaderItem(i,QTableWidgetItem(str(field+'('+str(i)+')')))
			if self.optimum['L'+field+'_OPT'][i] == 1:
				self.ui.TableOPTVarsProf.setItem(i,0,QTableWidgetItem('T'))
			else:
				self.ui.TableOPTVarsProf.setItem(i,0,QTableWidgetItem('F'))
			self.ui.TableOPTVarsProf.setItem(i,1,QTableWidgetItem(str(self.indata[field.lower()][i])))
			self.ui.TableOPTVarsProf.setItem(i,2,QTableWidgetItem(str(self.optimum['D'+field+'_OPT'][i])))
			self.ui.TableOPTVarsProf.setItem(i,3,QTableWidgetItem(str(self.optimum[field+'_MIN'][i])))
			self.ui.TableOPTVarsProf.setItem(i,4,QTableWidgetItem(str(self.optimum[field+'_MAX'][i])))
		self.ui.TableOPTVarsProf.blockSignals(False)

	def OPTVarsProf(self):
		# Handles changes in Profile table
		# Get the current profile
		data_name = self.ui.comboBoxStelVarsProfType.currentText()
		print(data_name)
		if data_name == 'Pressure (AM)':
			field = 'AM'
		elif data_name == 'Current (AC)':
			field = 'AC'
		elif data_name == 'Iota (AI)':
			field = 'AI'
		# Get changed item
		item = self.ui.TableOPTVarsProf.currentItem()
		if item is None:
			return
		col = self.ui.TableOPTVarsProf.currentColumn()
		row = self.ui.TableOPTVarsProf.currentRow()
		# Column determins value
		field = field.upper()
		if col == 0:
			field = 'L'+field+'_OPT'
			if item.text() == 'T':
				val = 1
			else:
				val = 0
		if col == 1:
			return
		elif col == 2:
			field = 'D'+field+'_OPT'
			val = float(item.text())
		elif col == 3:
			field = field+'_MIN'
			val = float(item.text())
		elif col == 4:
			field = field+'_MAX'
			val = float(item.text())
		print(field, val)
		temp = self.optimum[field]
		temp[row] = val
		self.optimum[field]=temp

	def UpdateOPTVarsExtcur(self):
		self.ui.TableOPTVarsExtcur.blockSignals(True)
		# See how many EXTCUR array vars there are
		if self.indata['lfreeb'] == 'T':
			nextcur = 0
			for i in range(99):
				if not (self.indata['extcur'] == 0):
					nextcur = i
			self.ui.TableOPTVarsExtcur.setRowCount(nextcur)
		else:
			self.ui.TableOPTVarsExtcur.setRowCount(1)
			return
		field = 'EXTCUR'
		for i in range(nextcur):
			if self.optimum['L'+field+'_OPT'][i] == 1:
				self.ui.TableOPTVarsProf.setItem(i,0,QTableWidgetItem('T'))
			else:
				self.ui.TableOPTVarsProf.setItem(i,0,QTableWidgetItem('F'))
			self.ui.TableOPTVarsProf.setItem(i,1,QTableWidgetItem(str(self.indata[field.lower()][i])))
			self.ui.TableOPTVarsProf.setItem(i,2,QTableWidgetItem(str(self.optimum['D'+field+'_OPT'][i])))
			self.ui.TableOPTVarsProf.setItem(i,3,QTableWidgetItem(str(self.optimum[field+'_MIN'][i])))
			self.ui.TableOPTVarsProf.setItem(i,4,QTableWidgetItem(str(self.optimum[field+'_MAX'][i])))
			self.ui.TableOPTVarsProf.setVerticalHeaderItem(i,QTableWidgetItem(str(field+'('+str(i)+')')))
		self.ui.TableOPTVarsProf.blockSignals(False)

	def OPTVarsExtcur(self):
		# Handles changes in Profile table
		# Get the current profile
		# Get changed item
		item = self.ui.TableOPTVarsExtcur.currentItem()
		if item is None:
			return
		col = self.ui.TableOPTVarsExtcur.currentColumn()
		row = self.ui.TableOPTVarsExtcur.currentRow()
		# Column determins value
		field = 'EXTCUR'
		if col == 0:
			field = 'L'+field+'_OPT'
			if item.text() == 'T':
				val = 1
			else:
				val = 0
		if col == 1:
			return
		elif col == 2:
			field = 'D'+field+'_OPT'
			val = float(item.text())
		elif col == 3:
			field = field+'_MIN'
			val = float(item.text())
		elif col == 4:
			field = field+'_MAX'
			val = float(item.text())
		print(field, val)
		temp = self.optimum[field]
		temp[row] = val
		self.optimum[field]=temp

	def UpdateOPTVarsBound(self):
		self.ui.TableOPTVarsBound.blockSignals(True)
		# Figure out which array to deal with
		data_name = self.ui.comboBoxStelVarsBoundType.currentText()
		if data_name == 'VMEC':
			field1 = 'BOUND'
			field2 = 'BOUND'
			field3 = 'BOUND'
			field4 = 'BOUND'
			nrows = (2*self.indata['ntor']+1)*(self.indata['mpol']-1)+self.indata['ntor']+1
			xm = np.ndarray((nrows,1))
			xn = np.ndarray((nrows,1))
			i=0
			for n in range(0,self.indata['ntor']+1):
				xm[i] = 0
				xn[i] = n
				i=i+1
			for m in range(1,self.indata['mpol']):
				for n in range(-self.indata['ntor'],self.indata['ntor']+1):
					if not(n<0 and m<1):
						xm[i] = m
						xn[i] = n
						i=i+1
		elif data_name == 'Hirshman-Breslau':
			field1 = 'RHO'
			field2 = 'RHO'
			field3 = 'BOUND'
			field4 = 'BOUND'
			nrows = (2*self.indata['ntor']+1)*(self.indata['mpol']-1)+self.indata['ntor']+1
			xm = np.ndarray((nrows,1))
			xn = np.ndarray((nrows,1))
			i=0
			for n in range(0,self.indata['ntor']+1):
				xm[i] = 0
				xn[i] = n
				i=i+1
			for m in range(1,self.indata['mpol']):
				for n in range(-self.indata['ntor'],self.indata['ntor']+1):
					if not(n<0 and m<1):
						xm[i] = m
						xn[i] = n
						i=i+1
		elif data_name == 'Garabedian':
			field1 = 'DELTAMN'
			field2 = 'DELTAMN'
			field3 = 'DELTA'
			field4 = 'DELTA'
			nrows = (2*self.indata['ntor']+1)*(2*self.indata['ntor']+1)
			xm = np.ndarray((nrows,1))
			xn = np.ndarray((nrows,1))
			i=0
			for m in range(-self.indata['mpol'],self.indata['mpol']+1):
				for n in range(-self.indata['ntor'],self.indata['ntor']+1):
						xm[i] = m
						xn[i] = n
						i=i+1
		self.ui.TableOPTVarsBound.setRowCount(nrows)
		for i in range(nrows):
			m = int(xm[i])
			n = int(xn[i])
			dex1 = n+101
			dex2 = m
			head_string = '('+str(n)+','+str(m)+')'
			self.ui.TableOPTVarsBound.setVerticalHeaderItem(i,QTableWidgetItem(head_string))
			if self.optimum['L'+field1+'_OPT'][dex1][dex2] == 1:
				self.ui.TableOPTVarsBound.setItem(i,0,QTableWidgetItem('T'))
			else:
				self.ui.TableOPTVarsBound.setItem(i,0,QTableWidgetItem('F'))
		#	self.ui.TableOPTVarsProf.setItem(i,1,QTableWidgetItem(str(self.indata[field.lower()][i])))
			self.ui.TableOPTVarsBound.setItem(i,2,QTableWidgetItem(str(self.optimum['D'+field2+'_OPT'][dex1][dex2])))
			self.ui.TableOPTVarsBound.setItem(i,3,QTableWidgetItem(str(self.optimum[field3+'_MIN'][dex1][dex2])))
			self.ui.TableOPTVarsBound.setItem(i,4,QTableWidgetItem(str(self.optimum[field4+'_MAX'][dex1][dex2])))
		self.ui.TableOPTVarsBound.blockSignals(False)

	def OPTVarsBound(self):
		# Handles changes in Boundary table
		# Get the current profile
		data_name = self.ui.comboBoxStelVarsBoundType.currentText()
		# Get changed item
		item = self.ui.TableOPTVarsBound.currentItem()
		if item is None:
			return
		col = self.ui.TableOPTVarsBound.currentColumn()
		row_text = self.ui.TableOPTVarsBound.verticalHeaderItem(col)
		row_text = row_text.text()
		row_text = row_text.replace('(','')
		row_text = row_text.replace(')','')
		n, m = row_text.split(',')
		n = int(n)
		m = int(m)
		if data_name == 'VMEC':
			field1 = 'BOUND'
			field2 = 'BOUND'
			field3 = 'BOUND'
			field4 = 'BOUND'
			dex1 = n+100
			dex2 = m

		# Column determins value
		if col == 0:
			field = 'L'+field1+'_OPT'
			if item.text() == 'T':
				val = 1
			else:
				val = 0
		if col == 1:
			return
		elif col == 2:
			field = 'D'+field2+'_OPT'
			val = float(item.text())
		elif col == 3:
			field = field3+'_MIN'
			val = float(item.text())
		elif col == 4:
			field = field4+'_MAX'
			val = float(item.text())
		print(field, val)
		temp = self.optimum[field]
		temp[dex1,dex2] = val
		self.optimum[field]=temp

	def LoadSTELLOPT(self):
		# Handles loading an stellopt file.
		w = QWidget()
		w.resize(320, 240)
		w.setWindowTitle("Hello World!")
		filename = QFileDialog.getOpenFileName(w, 'Open File', '.','STELLOPT (stellopt.*)')
		w.destroy
		# Read the file
		self.stel_data=read_stellopt(filename)
		self.optplot_list = ['ASPECT','BETA','CURTOR','EXTCUR','SEPARATRIX',\
					'PHIEDGE','RBTOR','R0','Z0','VOLUME','WP','KAPPA',\
					'B_PROBES','FARADAY','FLUXLOOPS','SEGROG','MSE',\
					'NE','NELINE','TE','TELINE','TI','TILINE',\
					'XICS','XICS_BRIGHT','XICS_W3','XICS_V','SXR','VPHI',\
					'IOTA','BALLOON','BOOTSTRAP','DKES','HELICITY','HELICITY_FULL',\
					'KINK','ORBIT','JDOTB','J_STAR','NEO','TXPORT','ECEREFLECT',\
					]
		self.ui.ComboBoxOPTplot_type.clear()
		self.ui.ComboBoxOPTplot_type.addItem('Chi-Squared')
		# Handle Chisquared plots
		for name in self.optplot_list:
			for item in self.stel_data:
				if (name+'_target' == item):
					self.ui.ComboBoxOPTplot_type.addItem(name)
		# Handle Special Plots
		self.ui.ComboBoxOPTplot_type.addItem('-----SPECIAL-----')
		for name in ['BALLOON','KINK','ORBIT','NEO','HELICITY','HELICITY_FULL',\
					'TXPORT','B_PROBES','FLUXLOOPS','SEGROG',\
					'NELINE','TELINE','TILINE',\
					'XICS','XICS_BRIGHT','XICS_W3','XICS_V',\
					'ECEREFLECT','SXR','IOTA','PRESS']:
			for item in self.stel_data:
				if (name+'_target' == item):
					self.ui.ComboBoxOPTplot_type.addItem(name+'_evolution')
		for name in ['NE','TE','TI','MSE']:
			for item in self.stel_data:
				if (name+'_target' == item):
					self.ui.ComboBoxOPTplot_type.addItem(name+'_evolution')
					self.ui.ComboBoxOPTplot_type.addItem(name+'_evolution_R')
					self.ui.ComboBoxOPTplot_type.addItem(name+'_evolution_Z')
		# Handle Wout Comparrison Plots
		self.workdir,ext = filename.split('stellopt.',1)
		files = os.listdir(self.workdir)
		if any('wout' in mystring for mystring in files):
			self.ui.ComboBoxOPTplot_type.addItem('----- VMEC -----')
			self.ui.ComboBoxOPTplot_type.addItem('Flux0')
			self.ui.ComboBoxOPTplot_type.addItem('FluxPI')
			self.ui.ComboBoxOPTplot_type.addItem('Pressure')
			self.ui.ComboBoxOPTplot_type.addItem('I-prime')
			self.ui.ComboBoxOPTplot_type.addItem('Current')
			self.ui.ComboBoxOPTplot_type.addItem('Iota')
			self.ui.ComboBoxOPTplot_type.addItem('q-prof')
			self.ui.ComboBoxOPTplot_type.addItem('<j*B>')
			wout_files = sorted([k for k in files if 'wout' in k])
			self.wout_files = sorted([k for k in wout_files if '_opt' not in k])
		# Handle Kinetic Profiles
		if any('tprof.' in mystring for mystring in files):
			self.ui.ComboBoxOPTplot_type.addItem('----- Kinetics -----')
			self.ui.ComboBoxOPTplot_type.addItem('Electron Temperature')
			self.ui.ComboBoxOPTplot_type.addItem('Electron Density')
			self.ui.ComboBoxOPTplot_type.addItem('Ion Temperature')
			self.ui.ComboBoxOPTplot_type.addItem('Z Effective')
			tprof_files = sorted([k for k in files if 'tprof.' in k])
			self.tprof_files = sorted([k for k in tprof_files if '_opt' not in k])
		# Handle Diagnostic Profiles
		if any('dprof.' in mystring for mystring in files):
			self.ui.ComboBoxOPTplot_type.addItem('----- Diagnostic -----')
			self.ui.ComboBoxOPTplot_type.addItem('XICS Emissivity')
			self.ui.ComboBoxOPTplot_type.addItem('E-Static Potential')
			self.dprof_files = sorted([k for k in files if 'dprof.' in k])
		# Handle Current Density Profiles
		if any('jprof.' in mystring for mystring in files):
			self.ui.ComboBoxOPTplot_type.addItem('----- Current Density -----')
			self.ui.ComboBoxOPTplot_type.addItem('Bootstrap Profile')
			self.ui.ComboBoxOPTplot_type.addItem('Beam Profile')
			self.ui.ComboBoxOPTplot_type.addItem('Total Current Profile')
			jprof_files = sorted([k for k in files if 'tprof.' in k])
			self.jprof_files = sorted([k for k in jprof_files if '_opt' not in k])
		

	def UpdateOptplot(self):
		# Handle plotting of 
		plot_name = self.ui.ComboBoxOPTplot_type.currentText()
		self.fig2.clf()
		#self.fig.delaxes(self.ax)
		#self.ax2 = self.fig2.add_subplot(111)
		self.ax2 = self.fig2.add_axes([0.2,0.2,0.7,0.7])
		if (plot_name == 'Chi-Squared'):
			chisq = ((self.stel_data['TARGETS'] - self.stel_data['VALS'])/self.stel_data['SIGMAS'])**2
			self.ax2.plot(self.stel_data['ITER'],np.sum(chisq,axis=1),'ok',label='Chisq Total')
			self.ax2.set_xlabel('Iteration')
			self.ax2.set_ylabel('Chi-Squared')
			self.ax2.set_title('Chi-Sqaured')
			self.ax2.set_yscale('log',basey=10)
			for name in self.optplot_list:
				if name+'_chisq' in self.stel_data:
					chisq_temp = self.stel_data[name+'_chisq']
					n = chisq_temp.shape;
					if (len(chisq_temp.shape) == 1):
						if n[0] > len(self.stel_data['ITER']):
							chisq_temp = np.sum(chisq_temp,axis=0)
					elif len(chisq_temp.shape) > 1:
						chisq_temp = np.sum(chisq_temp,axis=1)
					self.ax2.plot(self.stel_data['ITER'],chisq_temp,'o',fillstyle='none',label=name)
			self.ax2.legend()
		elif (plot_name in self.optplot_list):
			f = self.stel_data[plot_name+'_chisq']
			n = f.shape
			if (len(n)==0):
				# Single Time slice Single point
				n=0
			elif (len(n)==1):
				# Could be either mutli-time or single time
				if len(self.stel_data['ITER']) == n[0]:
					# Mutl-time single point
					n=0
				else:
					# Multi-channel single time
					f = np.sum(f,axis=0)
			else:
				# Multiple Time slices
				f = np.sum(f,axis=1)
			self.ax2.plot(self.stel_data['ITER'],f,'ok',fillstyle='none')
			self.ax2.set_xlabel('Iteration')
			self.ax2.set_ylabel('Chi-Squared')
			self.ax2.set_title(plot_name+' Chi-Sqaured')
			self.ax2.set_yscale('log',basey=10)
		elif (plot_name == 'BALLOON_evolution'):
			self.ax2.plot(self.stel_data['BALLOON_k'].T,self.stel_data['BALLOON_grate'].T,'o',fillstyle='none')
			self.ax2.set_xlabel('Radial Grid')
			self.ax2.set_ylabel('Growth Rate')
			self.ax2.set_title('COBRA Ballooning Stability (<0 Stable)')
		elif (plot_name == 'TXPORT_evolution'):
			self.ax2.plot(self.stel_data['TXPORT_s'].T,self.stel_data['TXPORT_equil'].T,'o',fillstyle='none')
			self.ax2.set_xlabel('Normalized Flux')
			self.ax2.set_ylabel('Proxy Function')
			self.ax2.set_title('Turbulent Transport Proxy')
			self.ax2.set_xlim((0,1))
		elif (plot_name == 'ORBIT_evolution'):
			self.ax2.plot(self.stel_data['ORBIT_s'].T,self.stel_data['ORBIT_equil'].T,'o',fillstyle='none')
			self.ax2.set_xlabel('Normalized Flux')
			self.ax2.set_ylabel('Orbit Losses')
			self.ax2.set_title('Gyro Particle Losses')
			self.ax2.set_xlim((0,1))
		elif (plot_name == 'NEO_evolution'):
			self.ax2.plot(self.stel_data['NEO_k'].T,self.stel_data['NEO_equil'].T,'o',fillstyle='none')
			self.ax2.set_xlabel('Radial Grid')
			self.ax2.set_ylabel('Epsilon Effective')
			self.ax2.set_title('Neoclassical Helical Ripple (NEO)')
		elif (plot_name == 'HELICITY_FULL_evolution'):
			self.ax2.plot(self.stel_data['HELICITY_FULL_equil'].T,'o',fillstyle='none')
			self.ax2.set_ylabel('Helicity')
			self.ax2.set_title('Boozer Spectrum Helicity')
		elif (plot_name == 'B_PROBE_evolution'):
			n=self.stel_data['B_PROBE_target'].shape
			x = np.ndarray((n[1],1))
			for j in range(n[1]): x[j]=j+1
			self.ax2.errorbar(x,self.stel_data['B_PROBE_target'].T,self.stel_data['B_PROBE_sigma'].T,marker='o',fillstyle='none')
			self.ax2.plot(x,self.stel_data['BPROBES_equil'].T,fillstyle='none')
			self.ax2.set_xlabel('B-Probe')
			self.ax2.set_ylabel('Signal')
			self.ax2.set_title('B-Probe Reconstruction')
		elif (plot_name == 'FLUXLOOPS_evolution'):
			n=self.stel_data['FLUXLOOPS_target'].shape
			y = self.stel_data['FLUXLOOPS_target'].T
			s = self.stel_data['FLUXLOOPS_sigma'].T
			e = self.stel_data['FLUXLOOPS_equil'].T
			b = s < 1E10
			dl = n[0]
			if len(n) > 1:
				b = b[:,0]
				x = np.ndarray((n[1],1))
				for j in range(len(x)): x[j]=j+1
				self.ax2.errorbar(x[b],y[b,0],s[b,0],fmt='sk',fillstyle='none')
				for l in range(dl):
					self.ax2.plot(x[b],e[b,l-1],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
			else:
				x = np.ndarray((n[0],1))
				for j in range(len(x)): x[j]=j+1
				self.ax2.errorbar(x[b],y[b],s[b],fmt='sk',fillstyle='none')
				self.ax2.plot(x[b],e[b],'o',fillstyle='none')
			self.ax2.set_xlabel('Fluxloop')
			self.ax2.set_ylabel('Signal')
			self.ax2.set_title('Fluxloop Reconstruction')
		elif (plot_name == 'SEGROG_evolution'):
			n=self.stel_data['SEGROG_target'].shape
			y = self.stel_data['SEGROG_target'].T
			s = self.stel_data['SEGROG_sigma'].T
			e = self.stel_data['SEGROG_equil'].T
			b = s < 1E10
			dl = n[0]
			if len(n) > 1:
				b = b[:,0]
				x = np.ndarray((n[1],1))
				for j in range(len(x)): x[j]=j+1
				self.ax2.errorbar(x[b],y[b,0],s[b,0],fmt='sk',fillstyle='none')
				for l in range(dl):
					self.ax2.plot(x[b],e[b,l-1],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
			else:
				x = np.ndarray((n[0],1))
				for j in range(len(x)): x[j]=j+1
				self.ax2.errorbar(x[b],y[b],s[b],fmt='sk',fillstyle='none')
				self.ax2.plot(x[b],e[b],'o',fillstyle='none')
			self.ax2.set_xlabel('Rogowski Coil')
			self.ax2.set_ylabel('Signal')
			self.ax2.set_title('Rogowski Reconstruction')
		elif (plot_name == 'ECEREFLECT_evolution'):
			n=self.stel_data['ECEREFLECT_target'].shape
			y = self.stel_data['ECEREFLECT_target'].T
			s = self.stel_data['ECEREFLECT_sigma'].T
			e = self.stel_data['ECEREFLECT_equil'].T
			dl = n[0]
			if (len(n)==0):
				# Single Time slice Single point
				x=np.ndarray((1,1))*0+1
				self.ax2.errorbar(x,y,s,fmt='sk',fillstyle='none')
				x = 1;
				self.ax2.plot(x,e,'o',fillstyle='none')
			elif (len(n)==1):
				# Could be either mutli-time or single time
				if len(self.stel_data['ITER']) == n[0]:
					# Mutl-time single point
					x = np.ndarray((n[0],1))*0+1
					self.ax2.errorbar(x[0],y[0],s[0],fmt='sk',fillstyle='none')
					self.ax2.plot(x,e,'o',fillstyle='none')
				else:
					# Multi-channel single time
					b = s < 1E10
					x = np.ndarray((n[0],1))
					for j in range(n[0]): x[j]=j+1
					self.ax2.errorbar(x[b],y[b],s[b],fmt='sk',fillstyle='none')
					self.ax2.plot(x[b],e[b],'o',fillstyle='none')
			else:
				# Multiple Time slices
				b = s[:,0] < 1E10
				x = np.ndarray((n[1],1))
				for j in range(n[1]): x[j]=j+1
				self.ax2.errorbar(x[b,0],y[b,0],s[b,0],fmt='sk',fillstyle='none')
				for l in range(dl):
					self.ax2.plot(x[b],e[b,l-1],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
			self.ax2.set_xlabel('ECE Channel')
			self.ax2.set_ylabel('Radiative Temp [eV]')
			self.ax2.set_title('ECE Reconstruction')
		elif (plot_name == 'KINK_evolution'):
			self.ax2.plot(self.stel_data['KINK_equil'],fmt='ok',fillstyle='none')
			self.ax2.plot(self.stel_data['KINK_target'],fmt='k')
			self.ax2.plot(self.stel_data['KINK_target']+self.stel_data['KINK_sigma'],fmt='k')
			self.ax2.plot(self.stel_data['KINK_target']-self.stel_data['KINK_sigma'],fmt='k')
			self.ax2.set_xlabel('Iteration')
			self.ax2.set_ylabel('???')
			self.ax2.set_title('?????KINK Evolution????')
		elif (plot_name == 'NE_evolution'):
			x = self.stel_data['NE_s'].T
			y = self.stel_data['NE_target'].T
			s = self.stel_data['NE_sigma'].T
			e = self.stel_data['NE_equil'].T
			n = y.shape
			if len(x.shape)>1:
				dl = n[1]
				self.ax2.errorbar(x[:,0],y[:,0],s[:,0],fmt='sk',fillstyle='none')
				for l in range(dl):
					self.ax2.plot(x[:,l-1],e[:,l-1],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
			else:
				self.ax2.errorbar(x[:],y[:],s[:],fmt='sk',fillstyle='none')
				self.ax2.plot(x[:],e[:],'o',fillstyle='none',color='g')
			self.ax2.set_xlabel('Normalized Flux')
			self.ax2.set_ylabel('Electron Density (norm)')
			self.ax2.set_title('Electron Density Reconstruction')
			self.ax2.set_xlim((0,1.6))
		elif (plot_name == 'NE_evolution_R'):
			x = self.stel_data['NE_R'].T
			y = self.stel_data['NE_target'].T
			s = self.stel_data['NE_sigma'].T
			e = self.stel_data['NE_equil'].T
			n = y.shape
			if len(x.shape)>1:
				dl = n[1]
				self.ax2.errorbar(x[:,0],y[:,0],s[:,0],fmt='sk',fillstyle='none')
				for l in range(dl):
					self.ax2.plot(x[:,l-1],e[:,l-1],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
			else:
				self.ax2.errorbar(x[:],y[:],s[:],fmt='sk',fillstyle='none')
				self.ax2.plot(x[:],e[:],'o',fillstyle='none',color='g')
			self.ax2.set_xlabel('R [m]')
			self.ax2.set_ylabel('Electron Density (norm)')
			self.ax2.set_title('Electron Density Reconstruction')
		elif (plot_name == 'NE_evolution_Z'):
			x = self.stel_data['NE_Z'].T
			y = self.stel_data['NE_target'].T
			s = self.stel_data['NE_sigma'].T
			e = self.stel_data['NE_equil'].T
			n = y.shape
			if len(x.shape)>1:
				dl = n[1]
				self.ax2.errorbar(x[:,0],y[:,0],s[:,0],fmt='sk',fillstyle='none')
				for l in range(dl):
					self.ax2.plot(x[:,l-1],e[:,l-1],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
			else:
				self.ax2.errorbar(x[:],y[:],s[:],fmt='sk',fillstyle='none')
				self.ax2.plot(x[:],e[:],'o',fillstyle='none',color='g')
			self.ax2.set_xlabel('Z [m]')
			self.ax2.set_ylabel('Electron Density (norm)')
			self.ax2.set_title('Electron Density Reconstruction')
		elif (plot_name == 'TE_evolution'):
			x=self.stel_data['TE_s'].T
			y=self.stel_data['TE_target'].T
			s=self.stel_data['TE_sigma'].T
			e = self.stel_data['TE_equil'].T
			n = y.shape
			if len(x.shape)>1:
				dl = n[1]
				self.ax2.errorbar(x[:,0],y[:,0],s[:,0],fmt='sk',fillstyle='none')
				for l in range(dl):
					self.ax2.plot(x[:,l-1],e[:,l-1],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
			else:
				self.ax2.errorbar(x[:],y[:],s[:],fmt='sk',fillstyle='none')
				self.ax2.plot(x[:],e[:],'o',fillstyle='none',color='g')
			self.ax2.set_xlabel('Normalized Flux')
			self.ax2.set_ylabel('Electron Temperature [keV]')
			self.ax2.set_title('Electron Temperature Reconstruction')
			self.ax2.set_xlim((0,1.6))
		elif (plot_name == 'TE_evolution_R'):
			x=self.stel_data['TE_R'].T
			y=self.stel_data['TE_target'].T
			s=self.stel_data['TE_sigma'].T
			e = self.stel_data['TE_equil'].T
			n = y.shape
			if len(x.shape)>1:
				dl = n[1]
				self.ax2.errorbar(x[:,0],y[:,0],s[:,0],fmt='sk',fillstyle='none')
				for l in range(dl):
					self.ax2.plot(x[:,l-1],e[:,l-1],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
			else:
				self.ax2.errorbar(x[:],y[:],s[:],fmt='sk',fillstyle='none')
				self.ax2.plot(x[:],e[:],'o',fillstyle='none',color='g')
			self.ax2.set_xlabel('R [m]')
			self.ax2.set_ylabel('Electron Temperature [keV]')
			self.ax2.set_title('Electron Temperature Reconstruction')
		elif (plot_name == 'TE_evolution_Z'):
			x=self.stel_data['TE_Z'].T
			y=self.stel_data['TE_target'].T
			s=self.stel_data['TE_sigma'].T
			e = self.stel_data['TE_equil'].T
			n = y.shape
			if len(x.shape)>1:
				dl = n[1]
				self.ax2.errorbar(x[:,0],y[:,0],s[:,0],fmt='sk',fillstyle='none')
				for l in range(dl):
					self.ax2.plot(x[:,l-1],e[:,l-1],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
			else:
				self.ax2.errorbar(x[:],y[:],s[:],fmt='sk',fillstyle='none')
				self.ax2.plot(x[:],e[:],'o',fillstyle='none',color='g')
			self.ax2.set_xlabel('Z [m]')
			self.ax2.set_ylabel('Electron Temperature [keV]')
			self.ax2.set_title('Electron Temperature Reconstruction')
		elif (plot_name == 'TI_evolution'):
			x=self.stel_data['TI_s'].T
			y=self.stel_data['TI_target'].T
			s=self.stel_data['TI_sigma'].T
			e = self.stel_data['TI_equil'].T
			n = y.shape
			if len(x.shape)>1:
				dl = n[1]
				self.ax2.errorbar(x[:,0],y[:,0],s[:,0],fmt='sk',fillstyle='none')
				for l in range(dl):
					self.ax2.plot(x[:,l-1],e[:,l-1],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
			else:
				self.ax2.errorbar(x[:],y[:],s[:],fmt='sk',fillstyle='none')
				self.ax2.plot(x[:],e[:],'o',fillstyle='none',color='g')
			self.ax2.set_xlabel('Normalized Flux')
			self.ax2.set_ylabel('Ion Temperature [keV]')
			self.ax2.set_title('Ion Temperature Reconstruction')
			self.ax2.set_xlim((0,1.6))
		elif (plot_name == 'TI_evolution_R'):
			x=self.stel_data['TI_R'].T
			y=self.stel_data['TI_target'].T
			s=self.stel_data['TI_sigma'].T
			e = self.stel_data['TI_equil'].T
			n = y.shape
			if len(x.shape)>1:
				dl = n[1]
				self.ax2.errorbar(x[:,0],y[:,0],s[:,0],fmt='sk',fillstyle='none')
				for l in range(dl):
					self.ax2.plot(x[:,l-1],e[:,l-1],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
			else:
				self.ax2.errorbar(x[:],y[:],s[:],fmt='sk',fillstyle='none')
				self.ax2.plot(x[:],e[:],'o',fillstyle='none',color='g')
			self.ax2.set_xlabel('R [m]')
			self.ax2.set_ylabel('Ion Temperature [keV]')
			self.ax2.set_title('Ion Temperature Reconstruction')
		elif (plot_name == 'TI_evolution_Z'):
			x=self.stel_data['TI_Z'].T
			y=self.stel_data['TI_target'].T
			s=self.stel_data['TI_sigma'].T
			e = self.stel_data['TI_equil'].T
			n = y.shape
			if len(x.shape)>1:
				dl = n[1]
				self.ax2.errorbar(x[:,0],y[:,0],s[:,0],fmt='sk',fillstyle='none')
				for l in range(dl):
					self.ax2.plot(x[:,l-1],e[:,l-1],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
			else:
				self.ax2.errorbar(x[:],y[:],s[:],fmt='sk',fillstyle='none')
				self.ax2.plot(x[:],e[:],'o',fillstyle='none',color='g')
			self.ax2.set_xlabel('Z [m]')
			self.ax2.set_ylabel('Ion Temperature [keV]')
			self.ax2.set_title('Ion Temperature Reconstruction')
		elif (plot_name == 'MSE_evolution'):
			x=self.stel_data['MSE_s'].T
			y=self.stel_data['MSE_target'].T
			s=self.stel_data['MSE_sigma'].T
			e = self.stel_data['MSE_equil'].T
			n = y.shape
			if len(x.shape)>1:
				dl = n[1]
				self.ax2.errorbar(x[:,0],y[:,0],s[:,0],fmt='sk',fillstyle='none')
				for l in range(dl):
					self.ax2.plot(x[:,l-1],e[:,l-1],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
			else:
				self.ax2.errorbar(x[:],y[:],s[:],fmt='sk',fillstyle='none')
				self.ax2.plot(x[:],e[:],'o',fillstyle='none',color='g')
			self.ax2.set_xlabel('Normalized Flux')
			self.ax2.set_ylabel('Pitch Angle')
			self.ax2.set_title('Motional Stark Effect')
			self.ax2.set_xlim((0,1.6))
		elif (plot_name == 'MSE_evolution_R'):
			x=self.stel_data['MSE_R'].T
			y=self.stel_data['MSE_target'].T
			s=self.stel_data['MSE_sigma'].T
			e = self.stel_data['MSE_equil'].T
			n = y.shape
			if len(x.shape)>1:
				dl = n[1]
				self.ax2.errorbar(x[:,0],y[:,0],s[:,0],fmt='sk',fillstyle='none')
				for l in range(dl):
					self.ax2.plot(x[:,l-1],e[:,l-1],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
			else:
				self.ax2.errorbar(x[:],y[:],s[:],fmt='sk',fillstyle='none')
				self.ax2.plot(x[:],e[:],'o',fillstyle='none',color='g')
			self.ax2.set_xlabel('R [m]')
			self.ax2.set_ylabel('Pitch Angle')
			self.ax2.set_title('Motional Stark Effect')
		elif (plot_name == 'MSE_evolution_Z'):
			x=self.stel_data['MSE_Z'].T
			y=self.stel_data['MSE_target'].T
			s=self.stel_data['MSE_sigma'].T
			e = self.stel_data['MSE_equil'].T
			n = y.shape
			if len(x.shape)>1:
				dl = n[1]
				self.ax2.errorbar(x[:,0],y[:,0],s[:,0],fmt='sk',fillstyle='none')
				for l in range(dl):
					self.ax2.plot(x[:,l-1],e[:,l-1],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
			else:
				self.ax2.errorbar(x[:],y[:],s[:],fmt='sk',fillstyle='none')
				self.ax2.plot(x[:],e[:],'o',fillstyle='none',color='g')
			self.ax2.set_xlabel('Z [m]')
			self.ax2.set_ylabel('Pitch Angle')
			self.ax2.set_title('Motional Stark Effect')
		elif (plot_name == 'IOTA_evolution'):
			x=self.stel_data['IOTA_s'].T
			y=self.stel_data['IOTA_target'].T
			s=self.stel_data['IOTA_sigma'].T
			e = self.stel_data['IOTA_equil'].T
			n = y.shape
			if len(x.shape)>1:
				dl = n[1]
				self.ax2.errorbar(x[:,0],y[:,0],s[:,0],fmt='sk',fillstyle='none')
				for l in range(dl):
					self.ax2.plot(x[:,l-1],e[:,l-1],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
			else:
				self.ax2.errorbar(x[:],y[:],s[:],fmt='sk',fillstyle='none')
				self.ax2.plot(x[:],e[:],'o',fillstyle='none',color='g')
			self.ax2.set_xlabel('Normalized Flux')
			self.ax2.set_ylabel('Iota')
			self.ax2.set_title('Rotational Transform')
		elif (plot_name == 'NELINE_evolution'):
			y=self.stel_data['NELINE_target'].T
			s=self.stel_data['NELINE_sigma'].T
			e = self.stel_data['NELINE_equil'].T
			n = y.shape
			if (len(n)==0):
				# Single Time slice Single point
				x=np.ndarray((1,1))*0+1
				self.ax2.errorbar(x,y,s,fmt='sk',fillstyle='none')
				self.ax2.plot(x,e,'og',fillstyle='none')
			elif (len(n)==1):
				# Could be either mutli-time or single time
				if len(self.stel_data['ITER']) == n[0]:
					# Mutl-time single point
					dl = n[0]
					x = np.ndarray((n[0],1))*0+1
					self.ax2.errorbar(x[0],y[0],s[0],fmt='sk',fillstyle='none')
					for l in range(dl): self.ax2.plot(x[l],e[l],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
				else:
					# Multi-channel single time
					x = np.ndarray((n[0],1))
					for j in range(n[0]): x[j]=j+1
					self.ax2.errorbar(x,y,s,fmt='sk',fillstyle='none')
					self.ax2.plot(x,e,'og',fillstyle='none')
			else:
				# Multiple Time slices
				dl = n[0]
				x = np.ndarray((n[1],1))
				for j in range(n[1]): x[j]=j+1
				self.ax2.errorbar(x[:,0],y[:,0],s[:,0],fmt='sk',fillstyle='none')
				for l in range(dl):
					self.ax2.plot(x,e[:,l-1],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
			self.ax2.set_xlabel('Channel')
			self.ax2.set_ylabel('Signal [m^{-2}]')
			self.ax2.set_title('Line-Int. Electron Density')
		elif (plot_name == 'TELINE_evolution'):
			y=self.stel_data['TELINE_target'].T
			s=self.stel_data['TELINE_sigma'].T
			e = self.stel_data['TELINE_equil'].T
			n = y.shape
			dl = n[0]
			if (len(n)==0):
				# Single Time slice Single point
				x=np.ndarray((1,1))*0+1
				self.ax2.errorbar(x,y,s,fmt='sk',fillstyle='none')
				self.ax2.plot(x,e,'og',fillstyle='none')
			elif (len(n)==1):
				# Could be either mutli-time or single time
				if len(self.stel_data['ITER']) == n[0]:
					# Mutl-time single point
					dl = n[0]
					x = np.ndarray((n[0],1))*0+1
					self.ax2.errorbar(x[0],y[0],s[0],fmt='sk',fillstyle='none')
					for l in range(dl): self.ax2.plot(x[l],e[l],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
				else:
					# Multi-channel single time
					x = np.ndarray((n[0],1))
					for j in range(n[0]): x[j]=j+1
					self.ax2.errorbar(x,y,s,fmt='sk',fillstyle='none')
					self.ax2.plot(x,e,'og',fillstyle='none')
			else:
				# Multiple Time slices
				dl = n[0]
				x = np.ndarray((n[1],1))
				for j in range(n[1]): x[j]=j+1
				self.ax2.errorbar(x[:,0],y[:,0],s[:,0],fmt='sk',fillstyle='none')
				for l in range(dl):
					self.ax2.plot(x,e[:,l-1],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
			self.ax2.set_xlabel('Channel')
			self.ax2.set_ylabel('Signal')
			self.ax2.set_title('Line-Int. Electron Temperature')
		elif (plot_name == 'TILINE_evolution'):
			y=self.stel_data['TILINE_target'].T
			s=self.stel_data['TILINE_sigma'].T
			e = self.stel_data['TILINE_equil'].T
			n = y.shape
			if (len(n)==0):
				# Single Time slice Single point
				x=np.ndarray((1,1))*0+1
				self.ax2.errorbar(x,y,s,fmt='sk',fillstyle='none')
				self.ax2.plot(x,e,'og',fillstyle='none')
			elif (len(n)==1):
				# Could be either mutli-time or single time
				if len(self.stel_data['ITER']) == n[0]:
					# Mutl-time single point
					dl = n[0]
					x = np.ndarray((n[0],1))*0+1
					self.ax2.errorbar(x[0],y[0],s[0],fmt='sk',fillstyle='none')
					for l in range(dl): self.ax2.plot(x[l],e[l],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
				else:
					# Multi-channel single time
					x = np.ndarray((n[0],1))
					for j in range(n[0]): x[j]=j+1
					self.ax2.errorbar(x,y,s,fmt='sk',fillstyle='none')
					self.ax2.plot(x,e,'og',fillstyle='none')
			else:
				# Multiple Time slices
				dl = n[0]
				x = np.ndarray((n[1],1))
				for j in range(n[1]): x[j]=j+1
				self.ax2.errorbar(x[:,0],y[:,0],s[:,0],fmt='sk',fillstyle='none')
				for l in range(dl):
					self.ax2.plot(x,e[:,l-1],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
			self.ax2.set_xlabel('Channel')
			self.ax2.set_ylabel('Signal')
			self.ax2.set_title('Line-Int. Ion Temperature')
		elif (plot_name == 'XICS_evolution'):
			y=self.stel_data['XICS_target'].T
			s=self.stel_data['XICS_sigma'].T
			e = self.stel_data['XICS_equil'].T
			n = y.shape
			if (len(n)==0):
				# Single Time slice Single point
				x=np.ndarray((1,1))*0+1
				self.ax2.errorbar(x,y,s,fmt='sk',fillstyle='none')
				self.ax2.plot(x,e,'og',fillstyle='none')
			elif (len(n)==1):
				# Could be either mutli-time or single time
				if len(self.stel_data['ITER']) == n[0]:
					# Mutl-time single point
					x = np.ndarray((n[0],1))*0+1
					self.ax2.errorbar(x[0],y[0],s[0],fmt='sk',fillstyle='none')
					for l in range(dl): self.ax2.plot(x[l],e[l],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
				else:
					# Multi-channel single time
					x = np.ndarray((n[0],1))
					for j in range(n[0]): x[j]=j+1
					self.ax2.errorbar(x,y,s,fmt='sk',fillstyle='none')
					self.ax2.plot(x,e,'og',fillstyle='none')
			else:
				# Multiple Time slices
				dl = n[1]
				x = np.ndarray((n[0],1))
				for j in range(n[0]): x[j]=j+1
				self.ax2.errorbar(x,y[:,0],s[:,0],fmt='sk',fillstyle='none')
				for l in range(dl):
					self.ax2.plot(x,e[:,l-1],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
			self.ax2.set_xlabel('Channel')
			self.ax2.set_ylabel('Signal [Arb.]')
			self.ax2.set_title('XICS Reconstruction')
		elif (plot_name == 'XICS_BRIGHT_evolution'):
			y=self.stel_data['XICS_BRIGHT_target'].T
			s=self.stel_data['XICS_BRIGHT_sigma'].T
			e = self.stel_data['XICS_BRIGHT_equil'].T
			n = y.shape
			if (len(n)==0):
				# Single Time slice Single point
				x=np.ndarray((1,1))*0+1
				self.ax2.errorbar(x,y,s,fmt='sk',fillstyle='none')
				self.ax2.plot(x,e,'og',fillstyle='none')
			elif (len(n)==1):
				# Could be either mutli-time or single time
				if len(self.stel_data['ITER']) == n[0]:
					# Mutl-time single point
					x = np.ndarray((n[0],1))*0+1
					self.ax2.errorbar(x[0],y[0],s[0],fmt='sk',fillstyle='none')
					for l in range(dl): self.ax2.plot(x[l],e[l],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
				else:
					# Multi-channel single time
					x = np.ndarray((n[0],1))
					for j in range(n[0]): x[j]=j+1
					self.ax2.errorbar(x,y,s,fmt='sk',fillstyle='none')
					self.ax2.plot(x,e,'og',fillstyle='none')
			else:
				# Multiple Time slices
				dl = n[1]
				x = np.ndarray((n[0],1))
				for j in range(n[0]): x[j]=j+1
				self.ax2.errorbar(x[:,0],y[:,0],s[:,0],fmt='sk',fillstyle='none')
				for l in range(dl):
					self.ax2.plot(x,e[:,l-1],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
			self.ax2.set_xlabel('Channel')
			self.ax2.set_ylabel('Signal [Arb.]')
			self.ax2.set_title('XICS Brightness Reconstruction')
		elif (plot_name == 'XICS_W3_evolution'):
			y=self.stel_data['XICS_W3_target'].T
			s=self.stel_data['XICS_W3_sigma'].T
			e = self.stel_data['XICS_W3_equil'].T
			n = y.shape
			if (len(n)==0):
				# Single Time slice Single point
				x=np.ndarray((1,1))*0+1
				self.ax2.errorbar(x,y,s,fmt='sk',fillstyle='none')
				self.ax2.plot(x,e,'og',fillstyle='none')
			elif (len(n)==1):
				# Could be either mutli-time or single time
				if len(self.stel_data['ITER']) == n[0]:
					# Mutl-time single point
					x = np.ndarray((n[0],1))*0+1
					self.ax2.errorbar(x[0],y[0],s[0],fmt='sk',fillstyle='none')
					for l in range(dl): self.ax2.plot(x[l],e[l],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
				else:
					# Multi-channel single time
					x = np.ndarray((n[0],1))
					for j in range(n[0]): x[j]=j+1
					self.ax2.errorbar(x,y,s,fmt='sk',fillstyle='none')
					self.ax2.plot(x,e,'og',fillstyle='none')
			else:
				# Multiple Time slices
				dl = n[1]
				x = np.ndarray((n[0],1))
				for j in range(n[0]): x[j]=j+1
				self.ax2.errorbar(x[:,0],y[:,0],s[:,0],fmt='sk',fillstyle='none')
				for l in range(dl):
					self.ax2.plot(x,e[:,l-1],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
			self.ax2.set_xlabel('Channel')
			self.ax2.set_ylabel('Signal [Arb.]')
			self.ax2.set_title('XICS W3 (Te) Reconstruction')
		elif (plot_name == 'XICS_V_evolution'):
			y=self.stel_data['XICS_V_target'].T
			s=self.stel_data['XICS_V_sigma'].T
			e = self.stel_data['XICS_V_equil'].T
			n = y.shape
			if (len(n)==0):
				# Single Time slice Single point
				x=np.ndarray((1,1))*0+1
				self.ax2.errorbar(x,y,s,fmt='sk',fillstyle='none')
				self.ax2.plot(x,e,'og',fillstyle='none')
			elif (len(n)==1):
				# Could be either mutli-time or single time
				if len(self.stel_data['ITER']) == n[0]:
					# Mutl-time single point
					x = np.ndarray((n[0],1))*0+1
					self.ax2.errorbar(x[0],y[0],s[0],fmt='sk',fillstyle='none')
					for l in range(dl): self.ax2.plot(x[l],e[l],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
				else:
					# Multi-channel single time
					x = np.ndarray((n[0],1))
					for j in range(n[0]): x[j]=j+1
					self.ax2.errorbar(x,y,s,fmt='sk',fillstyle='none')
					self.ax2.plot(x,e,'og',fillstyle='none')
			else:
				# Multiple Time slices
				dl = n[1]
				x = np.ndarray((n[0],1))
				for j in range(n[0]): x[j]=j+1
				self.ax2.errorbar(x[:,0],y[:,0],s[:,0],fmt='sk',fillstyle='none')
				for l in range(dl):
					self.ax2.plot(x,e[:,l-1],'o',fillstyle='none',color=_plt.cm.brg(l/(dl-1)))
			self.ax2.set_xlabel('Channel')
			self.ax2.set_ylabel('Signal [Arb.]')
			self.ax2.set_title('XICS Velocity Reconstruction')
		elif (plot_name == 'Pressure'):
			l=0
			dl = len(self.wout_files)-1
			if dl == 0 : dl = 1 
			for string in self.wout_files:
				if 'wout' in string:
					vmec_data=read_vmec(self.workdir+string)
					ns = vmec_data['ns']
					nflux = np.ndarray((ns,1))
					for j in range(ns): nflux[j]=j/(ns-1)
					self.ax2.plot(nflux,vmec_data['presf']/1000,color=_plt.cm.brg(l/dl))
					l=l+1
			self.ax2.set_xlabel('Norm Tor. Flux (s)')
			self.ax2.set_ylabel('Pressure [kPa]')
			self.ax2.set_title('VMEC Pressure Evolution')
			self.ax2.set_xlim((0,1))
		elif (plot_name == 'I-prime'):
			l=0
			dl = len(self.wout_files)-1
			if dl == 0 : dl = 1 
			for string in self.wout_files:
				if 'wout' in string:
					vmec_data=read_vmec(self.workdir+string)
					ns = vmec_data['ns']
					nflux = np.ndarray((ns,1))
					for j in range(ns): nflux[j]=j/(ns-1)
					self.ax2.plot(nflux,vmec_data['jcurv']/1000,color=_plt.cm.brg(l/dl))
					l=l+1
			self.ax2.set_xlabel('Norm Tor. Flux (s)')
			self.ax2.set_ylabel('Current Density [kA/m^{-2}]')
			self.ax2.set_title('VMEC Current Density Evolution')
			self.ax2.set_xlim((0,1))
		elif (plot_name == 'Iota'):
			l=0
			dl = len(self.wout_files)-1
			for string in self.wout_files:
				if 'wout' in string:
					vmec_data=read_vmec(self.workdir+string)
					ns = vmec_data['ns']
					nflux = np.ndarray((ns,1))
					for j in range(ns): nflux[j]=j/(ns-1)
					self.ax2.plot(nflux,vmec_data['iotaf'],color=_plt.cm.brg(l/dl))
					l=l+1
			self.ax2.set_xlabel('Norm Tor. Flux (s)')
			self.ax2.set_ylabel('Rotational Transform \iota')
			self.ax2.set_title('VMEC Rotational Transform Evolution')
			self.ax2.set_xlim((0,1))
		elif (plot_name == 'q-prof'):
			l=0
			dl = len(self.wout_files)
			for string in self.wout_files:
				if 'wout' in string:
					vmec_data=read_vmec(self.workdir+string)
					ns = vmec_data['ns']
					nflux = np.ndarray((ns,1))
					for j in range(ns): nflux[j]=j/(ns-1)
					self.ax2.plot(nflux,1.0/vmec_data['iotaf'],color=_plt.cm.brg(l/dl))
					l=l+1
			self.ax2.set_xlabel('Norm Tor. Flux (s)')
			self.ax2.set_ylabel('Safety Factor q')
			self.ax2.set_title('VMEC Safety Factor Evolution')
			self.ax2.set_xlim((0,1))
		elif (plot_name == 'Current'):
			l=0
			dl = len(self.wout_files)-1
			for string in self.wout_files:
				if 'wout' in string:
					vmec_data=read_vmec(self.workdir+string)
					ns = vmec_data['ns']
					nflux = np.ndarray((ns,1))
					for j in range(ns): nflux[j]=j/(ns-1)
					self.ax2.plot(nflux,np.cumsum(vmec_data['jcurv']),color=_plt.cm.brg(l/dl))
					l=l+1
			self.ax2.set_xlabel('Norm Tor. Flux (s)')
			self.ax2.set_ylabel('Current [kA]')
			self.ax2.set_title('VMEC Current Evolution')
			self.ax2.set_xlim((0,1))
		elif (plot_name == '<j*B>'):
			l=0
			dl = len(self.wout_files)-1
			for string in self.wout_files:
				if 'wout' in string:
					vmec_data=read_vmec(self.workdir+string)
					ns = vmec_data['ns']
					nflux = np.ndarray((ns,1))
					for j in range(ns): nflux[j]=j/(ns-1)
					self.ax2.plot(nflux,vmec_data['jdotb'],color=_plt.cm.brg(l/dl))
					l=l+1
			self.ax2.set_xlabel('Norm Tor. Flux (s)')
			self.ax2.set_ylabel('<j.B>')
			self.ax2.set_title('VMEC <j.B> Evolution')
			self.ax2.set_xlim((0,1))
		elif (plot_name == 'Flux0'):
			l=0
			dl = len(self.wout_files)-1
			for string in self.wout_files:
				if 'wout' in string:
					vmec_data=read_vmec(self.workdir+string)
					ns = vmec_data['ns']
					nu = 256
					nv = 1
					nflux = np.ndarray((ns,1))
					theta = np.ndarray((nu,1))
					zeta = np.ndarray((nv,1))
					for j in range(ns): nflux[j]=j/(ns-1)
					for j in range(nu): theta[j]=2*pi*j/(nu-1)
					zeta[0]=0;
					r=cfunct(theta,zeta,vmec_data['rmnc'],vmec_data['xm'],vmec_data['xn'])
					z=sfunct(theta,zeta,vmec_data['zmns'],vmec_data['xm'],vmec_data['xn'])
					self.ax2.plot(r[ns-1,:,0],z[ns-1,:,0],color=_plt.cm.brg(l/dl))
					self.ax2.plot(r[0,0,0],z[0,0,0],'+',color=_plt.cm.brg(l/dl))
					l=l+1
			self.ax2.set_xlabel('R [m]')
			self.ax2.set_ylabel('Z [m]')
			self.ax2.set_title('VMEC Flux Surface Evolution (phi=0)')
			self.ax2.set_aspect('equal')
		elif (plot_name == 'FluxPI'):
			l=0
			dl = len(self.wout_files)-1
			if dl == 0 : dl = 1 
			for string in self.wout_files:
				if 'wout' in string:
					vmec_data=read_vmec(self.workdir+string)
					ns = vmec_data['ns']
					nu = 256
					nv = 1
					nflux = np.ndarray((ns,1))
					theta = np.ndarray((nu,1))
					zeta = np.ndarray((nv,1))
					for j in range(ns): nflux[j]=j/(ns-1)
					for j in range(nu): theta[j]=2*pi*j/(nu-1)
					zeta[0]=pi/vmec_data['nfp'];
					r=cfunct(theta,zeta,vmec_data['rmnc'],vmec_data['xm'],vmec_data['xn'])
					z=sfunct(theta,zeta,vmec_data['zmns'],vmec_data['xm'],vmec_data['xn'])
					self.ax2.plot(r[ns-1,:,0],z[ns-1,:,0],color=_plt.cm.brg(l/dl))
					self.ax2.plot(r[0,0,0],z[0,0,0],'+',color=_plt.cm.brg(l/dl))
					l=l+1
			self.ax2.set_xlabel('R [m]')
			self.ax2.set_ylabel('Z [m]')
			self.ax2.set_title('VMEC Flux Surface Evolution (phi=0)')
			self.ax2.set_aspect('equal')
		elif (plot_name == 'Electron Temperature'):
			l=0
			dl = len(self.tprof_files)-1
			if dl == 0 : dl = 1 
			for string in self.tprof_files:
				if 'tprof' in string:
					tprof = np.loadtxt(self.workdir+string,skiprows=1)
					self.ax2.plot(tprof[:,0],tprof[:,2]/1E3,color=_plt.cm.brg(l/dl))
					l=l+1
			self.ax2.set_xlabel('Norm Tor. Flux (s)')
			self.ax2.set_ylabel('Temperature [keV]')
			self.ax2.set_title('Electron Temperature Profile')
			self.ax2.set_xlim((0,1))
		elif (plot_name == 'Electron Density'):
			l=0
			dl = len(self.tprof_files)-1
			if dl == 0 : dl = 1 
			for string in self.tprof_files:
				if 'tprof' in string:
					tprof = np.loadtxt(self.workdir+string,skiprows=1)
					self.ax2.plot(tprof[:,0],tprof[:,1]/1E19,color=_plt.cm.brg(l/dl))
					l=l+1
			self.ax2.set_xlabel('Norm Tor. Flux (s)')
			self.ax2.set_ylabel('Density [m^{-3}]')
			self.ax2.set_title('Electron Density Profile')
			self.ax2.set_xlim((0,1))
		elif (plot_name == 'Ion Temperature'):
			l=0
			dl = len(self.tprof_files)-1
			if dl == 0 : dl = 1 
			for string in self.tprof_files:
				if 'tprof' in string:
					tprof = np.loadtxt(self.workdir+string,skiprows=1)
					self.ax2.plot(tprof[:,0],tprof[:,3]/1E3,color=_plt.cm.brg(l/dl))
					l=l+1
			self.ax2.set_xlabel('Norm Tor. Flux (s)')
			self.ax2.set_ylabel('Temperature [keV]')
			self.ax2.set_title('Ion Temperature Profile')
			self.ax2.set_xlim((0,1))
		elif (plot_name == 'Z Effective'):
			l=0
			dl = len(self.tprof_files)-1
			if dl == 0 : dl = 1 
			for string in self.tprof_files:
				if 'tprof' in string:
					tprof = np.loadtxt(self.workdir+string,skiprows=1)
					self.ax2.plot(tprof[:,0],tprof[:,4],color=_plt.cm.brg(l/dl))
					l=l+1
			self.ax2.set_xlabel('Norm Tor. Flux (s)')
			self.ax2.set_ylabel('Z_{Eff}')
			self.ax2.set_title('Z Effective Profile')
			self.ax2.set_xlim((0,1))
		elif (plot_name == 'XICS Emissivity'):
			l=0
			dl = len(self.dprof_files)-1
			if dl == 0 : dl = 1 
			for string in self.dprof_files:
				if 'dprof' in string:
					dprof = np.loadtxt(self.workdir+string,skiprows=1)
					self.ax2.plot(dprof[:,0],dprof[:,1],color=_plt.cm.brg(l/dl))
					l=l+1
			self.ax2.set_xlabel('Norm Tor. Flux (s)')
			self.ax2.set_ylabel('Effective Emissivity')
			self.ax2.set_title('XICS Emissivity Profile')
			self.ax2.set_xlim((0,1))
		elif (plot_name == 'E-Static Potential'):
			l=0
			dl = len(self.dprof_files)
			for string in self.dprof_files:
				if 'dprof' in string:
					dprof = np.loadtxt(self.workdir+string,skiprows=1)
					self.ax2.plot(dprof[:,0],dprof[:,2]/1000,color=_plt.cm.brg(l/dl))
					l=l+1
			self.ax2.set_xlabel('Norm Tor. Flux (s)')
			self.ax2.set_ylabel('E-Static Potential (kV)')
			self.ax2.set_title('Electrostatic Potential Profile')
			self.ax2.set_xlim((0,1))
		elif (plot_name == 'Bootstrap Profile'):
			l=0
			dl = len(self.jprof_files)-1
			if dl == 0 : dl = 1 
			for string in self.jprof_files:
				if 'jprof' in string:
					jprof = np.loadtxt(self.workdir+string,skiprows=1)
					self.ax2.plot(jprof[:,0],jprof[:,2],color=_plt.cm.brg(l/dl))
					l=l+1
			self.ax2.set_xlabel('Norm Tor. Flux (s)')
			self.ax2.set_ylabel('Current Density [kA/m^-2]')
			self.ax2.set_title('Bootstrap Current Profile')
			self.ax2.set_xlim((0,1))
		elif (plot_name == 'Beam Profile'):
			l=0
			dl = len(self.jprof_files)-1
			if dl == 0 : dl = 1 
			for string in self.jprof_files:
				if 'jprof' in string:
					jprof = np.loadtxt(self.workdir+string,skiprows=1)
					self.ax2.plot(jprof[:,0],jprof[:,1],color=_plt.cm.brg(l/dl))
					l=l+1
			self.ax2.set_xlabel('Norm Tor. Flux (s)')
			self.ax2.set_ylabel('Current Density [kA/m^-2]')
			self.ax2.set_title('Beam Current Profile')
			self.ax2.set_xlim((0,1))
		elif (plot_name == 'Total Current Profile'):
			l=0
			dl = len(self.jprof_files)-1
			if dl == 0 : dl = 1 
			jprof = np.loadtxt(self.workdir+self.jprof_files[0],skiprows=1)
			self.ax2.plot(jprof[:,0],jprof[:,2],'--b')
			self.ax2.plot(jprof[:,0],jprof[:,1],':b')
			self.ax2.plot(jprof[:,0],jprof[:,3],'b')
			jprof = np.loadtxt(self.workdir+self.jprof_files[dl-1],skiprows=1)
			self.ax2.plot(jprof[:,0],jprof[:,2],'--g')
			self.ax2.plot(jprof[:,0],jprof[:,1],':g')
			self.ax2.plot(jprof[:,0],jprof[:,3],'g')
			self.ax2.set_xlabel('Norm Tor. Flux (s)')
			self.ax2.set_ylabel('Current Density [kA/m^-2]')
			self.ax2.set_title('Total Current Profile')
			self.ax2.set_xlim((0,1))
		self.canvas2.draw()




if __name__ == "__main__":
	app = QApplication(sys.argv) 
	window = MyApp() 
	window.show() 
	sys.exit(app.exec_())
