#!/usr/bin/env python3
import sys, os
os.environ['ETS_TOOLKIT'] = 'qt5'
import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as _plt
import numpy as np                    #For Arrays
from math import pi
#QT5
from PyQt5 import uic, QtGui, QtWidgets
from PyQt5.QtWidgets import QMainWindow, QApplication, QVBoxLayout, QSizePolicy, QWidget, QFileDialog, QTableWidgetItem
from PyQt5.QtGui import QIcon
# Matplotlib
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits import mplot3d
# VTK
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
#
from libstell import vmec
from libstell import stellopt
from libstell import plot3D

try:
	qtCreatorPath=os.environ["STELLOPT_PATH"]
except KeyError:
	print("Please set environment variable STELLOPT_PATH")
	sys.exit(1)

qtCreatorFile = os.path.join(qtCreatorPath,'pySTEL','STELLOPT.ui')
#qtCreatorFile = qtCreatorPath+"/pySTEL/STELLOPT.ui" # Enter file here.
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

class MyApp(QMainWindow):
	def __init__(self):
		super(MyApp, self).__init__()
		self.ui = Ui_MainWindow()
		self.ui.setupUi(self) 
		self.setStyleSheet("background-color: white;")
		# Init
		self.indata    = vmec.VMEC_INDATA()
		self.optimum   = stellopt.STELLOPT_INPUT()
		self.stel_data = stellopt.STELLOPT()
		# Set default values
		self.indata.read_indata('')
		#self.optimum.read_input('nofile')
		self.ui.tabMain.setTabEnabled(0,False)
		#self.ui.tabMain.setTabEnabled(2,False)
		#self.ui.tabStelVars.setTabEnabled(3,False)


		# Setup Defaults
		#iunit = 55
		#istat = 0
		#self.indata=read_indata_namelist(iunit,istat) # dummy just to get going
		#self.indata['ntor']=4
		#self.indata['mpol']=9
		#self.indata['nfp']=3
		#self.ui.tabMain.setTabEnabled(0,False)
		#self.ui.tabMain.setTabEnabled(2,False)
		# Set the OPTIMUM DEFAULTS (will repalce with something like read_indata_namelist)
		#self.optimum=read_stellopt_namelist(iunit,istat)
		# Setup Components
		#self.ui.TableArrays.setRowCount(1)
		#self.ui.TableArrays.setColumnCount(20)
		#for num,item in enumerate(self.indata['am'], start=0):
		#	self.ui.TableArrays.setItem(0,num, QTableWidgetItem(str(item)))
		# Setup Plot VMEC
		self.fig = Figure(figsize=(2,2),dpi=100)
		self.ax = self.fig.add_subplot(111)
		self.canvas = FigureCanvas(self.fig)
		self.ui.Plotbox.addWidget(self.canvas)
		# VTK stuff
		self.frame_vtk = QWidget()
		self.vtkWidget = QVTKRenderWindowInteractor(self.frame_vtk)
		self.ui.Plotbox.addWidget(self.vtkWidget)
		self.vtkWidget.Initialize()
		self.vtkWidget.Start()
		# Create a VTK renderer and add it to the render window
		self.plt_vmec = plot3D.PLOT3D(lwindow=False)
		self.vtkWidget.GetRenderWindow().AddRenderer(self.plt_vmec.renderer)
		self.vtkWidget.hide()
		# Setup Plot STELLOPT
		self.fig2 = Figure(figsize=(2,2),dpi=100)
		self.ax2 = self.fig2.add_subplot(111)
		self.canvas2 = FigureCanvas(self.fig2)
		self.ui.OPTplot_box.addWidget(self.canvas2)
		self.toolbar = NavigationToolbar(self.canvas2, self)
		self.ui.OPTplot_box.addWidget(self.toolbar)
		# Setup STELLOPT Pannels
		#self.UpdateOPTtype()
		#self.UpdateOPTVarsScalar()
		#self.UpdateOPTVarsProf()
		#self.UpdateOPTVarsExtcur()
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
		self.ui.ButtonLoadOptimum.clicked.connect(self.LoadOptimum)
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
		self.ui.ButtonPlotSTELLOPT.clicked.connect(self.PlotSTELLOPT)

	def UpdateMpol(self):
		strtmp = self.ui.TextMpol.text()
		inttmp = int(strtmp)
		setattr(self.indata,'mpol',inttmp)
		#self.indata['mpol'] = inttmp
		#set_module_var('vmec_input','mpol',inttmp)
		return

	def UpdateNtor(self):
		strtmp = self.ui.TextNtor.text()
		inttmp = int(strtmp)
		setattr(self.indata,'ntor',inttmp)
		#self.indata['ntor'] = inttmp
		#set_module_var('vmec_input','ntor',inttmp)
		return
		return

	def UpdateNfp(self):
		strtmp = self.ui.TextNfp.text()
		inttmp = int(strtmp)
		setattr(self.indata,'nfp',inttmp)
		#self.indata['nfp'] = inttmp
		#set_module_var('vmec_input','nfp',inttmp)
		return

	def LoadIndata(self):
		# Handles loading an indata file.
		w = QWidget()
		w.resize(320, 240)
		w.setWindowTitle("Load VMEC INDATA")
		filename = QFileDialog.getOpenFileName(w, 'Open File', '.')
		w.destroy
		# Now read the file
		iunit = 27
		istat = 0
		recl  = 1
		#temp=safe_open(iunit,istat,filename,'old','formatted',recl,'sequential','none')
		self.indata.read_indata(filename[0])
		#self.indata=read_indata_namelist(iunit,istat)
		#safe_close(iunit)
		# Now update the UI
		self.ui.TextTcon0.setText(str(self.indata.tcon0)) 
		self.ui.TextDelt.setText(str(self.indata.delt)) 
		self.ui.TextMpol.setText(str(self.indata.mpol)) 
		self.ui.TextNtor.setText(str(self.indata.ntor)) 
		self.ui.TextNtheta.setText(str(self.indata.ntheta)) 
		self.ui.TextNzeta.setText(str(self.indata.nzeta)) 
		self.ui.TextNfp.setText(str(self.indata.nfp)) 
		self.ui.TextNvacSkip.setText(str(self.indata.nvacskip)) 
		self.ui.TextMgridFile.setText(self.indata.mgrid_file) 
		self.ui.CheckboxLfreeb.setChecked(self.indata.lfreeb)
		self.ui.TextGamma.setText(str(self.indata.gamma))
		self.ui.TextPhiEdge.setText(str(self.indata.phiedge))
		self.ui.TextCurtor.setText(str(self.indata.curtor))
		self.ui.TextPresScale.setText(str(self.indata.pres_scale)) 
		self.ui.TextBloat.setText(str(self.indata.bloat)) 
		self.ui.TextSpresPed.setText(str(self.indata.spres_ped))
		self.ui.ComboBoxNcurr.setCurrentIndex(self.indata.ncurr)
		# Handle the ComboBoxPType
		# Handle NS Array table widget
		ns_mask=self.indata.ns_array!=0
		ns_array = self.indata.ns_array[ns_mask]
		ftol_array = self.indata.ftol_array[ns_mask]
		niter_array = self.indata.niter_array[ns_mask]
		self.ui.TableNsArray.setColumnCount(len(ns_array))
		for num,item in enumerate(ns_array, start=0):
			self.ui.TableNsArray.setItem(0,num, QTableWidgetItem(str(item)))
			self.ui.TableNsArray.setItem(1,num, QTableWidgetItem(str(niter_array[num])))
			self.ui.TableNsArray.setItem(2,num, QTableWidgetItem(str(ftol_array[num])))
		self.ui.TableNsArray.show()
		self.UpdateArrays()
		#self.ui.tabStelVars.setTabEnabled(3,True)

	def WriteIndata(self):
		# Handles loading an indata file.
		w = QWidget()
		w.resize(320, 240)
		w.setWindowTitle("Save VMEC INDATA")
		filename = QFileDialog.getSaveFileName(w, 'Open File', '.')
		w.destroy
		# Update the module with all the not-updated values
		setattr(self.indata,'tcon0',float(self.ui.TextTcon0.text()))
		setattr(self.indata,'delt',float(self.ui.TextDelt.text()))
		setattr(self.indata,'ntheta',int(self.ui.TextNtheta.text()))
		setattr(self.indata,'nzeta',int(self.ui.TextNzeta.text()))
		setattr(self.indata,'nfp',int(self.ui.TextNfp.text()))
		setattr(self.indata,'nvacskip',int(self.ui.TextNvacSkip.text()))
		setattr(self.indata,'gamma',float(self.ui.TextGamma.text()))
		setattr(self.indata,'phiedge',float(self.ui.TextPhiEdge.text()))
		setattr(self.indata,'pres_scale',float(self.ui.TextPresScale.text()))
		setattr(self.indata,'bloat',float(self.ui.TextBloat.text()))
		setattr(self.indata,'spres_ped',float(self.ui.TextSpresPed.text()))
		setattr(self.indata,'curtor',float(self.ui.TextCurtor.text()))
		setattr(self.indata,'mgrid_file',self.ui.TextMgridFile.text())
		# Update NCURR
		dex = self.ui.ComboBoxNcurr.currentIndex()
		setattr(self.indata,'ncurr',int(dex))
		# Now write the file
		self.indata.write_indata(filename[0])

	def UpdatePType(self):
		data_name=self.ui.ComboBoxArrays.currentText()
		data_name = data_name.lower()
		type_name=self.ui.ComboBoxPType.currentText()
		type_name = type_name.lower()
		if data_name == 'am' or data_name == 'am_aux':
			setattr(self.indata,'pmass_type',type_name)
		elif data_name == 'ac' or data_name == 'ac_aux':
			setattr(self.indata,'pcurr_type',type_name)
		elif data_name == 'ai' or data_name == 'ai_aux':
			setattr(self.indata,'piota_type',type_name)
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
			temp = getattr(self.indata,'ns_array')
			temp[col] = int(item.text())
			setattr(self.indata,'ns_array',temp)
		elif row == 1: #NITER_ARRAY
			temp = getattr(self.indata,'niter_array')
			temp[col] = int(item.text())
			setattr(self.indata,'niter_array',temp)
		elif row ==2: #FTOL
			temp = getattr(self.indata,'ftol_array')
			temp[col] = float(item.text())
			setattr(self.indata,'ftol_array',temp)

	def UpdateArrays(self):
		# Updates values in array box
		dex=self.ui.ComboBoxArrays.currentIndex()
		data_name=self.ui.ComboBoxArrays.currentText()
		data_name = data_name.lower()
		# Update PType accordingly
		strtmp = 'none'
		if data_name == 'am':
			strtmp = str(self.indata.pmass_type).strip()
		if data_name == 'ai':
			strtmp = str(self.indata.piota_type).strip()
		if data_name == 'ac' or data_name == 'ac_aux':
			strtmp = str(self.indata.pcurr_type).strip()
		if data_name == 'am_aux':
			setattr(self.indata,'pmass_type','akima_spline')
			strtmp = 'akima_spline'
		if data_name == 'ai_aux':
			setattr(self.indata,'piota_type','akima_spline')
			strtmp = 'akima_spline'
		if strtmp != 'none':
			self.ui.ComboBoxPType.setCurrentIndex(self.ui.ComboBoxPType.findText(strtmp))
		# Update the table
		self.ui.TableArrays.setRowCount(1)
		self.ui.TableArrays.setColumnCount(20)
		#self.ui.TableArrays.setVerticalHeaderLabels('1')
		self.ui.TableArrays.setHorizontalHeaderLabels('0;1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20'.split(';'))
		if data_name in ['am','ac','ai']:
			temp = getattr(self.indata,data_name)
			for num,item in enumerate(temp, start=0):
				self.ui.TableArrays.setItem(0,num, QTableWidgetItem(str(item)))
		elif data_name in ['am_aux','ac_aux','ai_aux']:
			temp = getattr(self.indata,data_name+'_s')
			aux_mask = temp>=0
			self.ui.TableArrays.setRowCount(2)
			self.ui.TableArrays.setColumnCount(50)
			if (np.count_nonzero(aux_mask) > 0):
				temp_val = getattr(self.indata,data_name+'_f')
				self.ui.TableArrays.setColumnCount(np.count_nonzero(aux_mask))
				for num,item in enumerate(temp[aux_mask], start=0):
					self.ui.TableArrays.setItem(0,num, QTableWidgetItem(str(item)))
					self.ui.TableArrays.setItem(1,num, QTableWidgetItem(str(temp_val[num])))
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
				setattr(self.indata,data_name+'_s',np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]))
				setattr(self.indata,data_name+'_f',np.array([1.0, 0.8, 0.6, 0.4, 0.2, 0.0]))
				#self.indata[data_name+'_s']=np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
				#self.indata[data_name+'_f']=np.array([1.0, 0.8, 0.6, 0.4, 0.2, 0.0])
				#set_module_var('vmec_input',data_name+'_s',self.indata[data_name+'_s'])
				#set_module_var('vmec_input',data_name+'_f',self.indata[data_name+'_f'])
		elif data_name == 'rbc' or data_name == 'zbs' or data_name == 'rbs' or data_name == 'zbc':
			self.ui.TableArrays.blockSignals(True)
			# Set Array Size
			self.ui.TableArrays.setRowCount(self.indata.ntor*2+1)
			self.ui.TableArrays.setColumnCount(self.indata.mpol)
			# Set Headings
			temp =''
			for i in range(2*self.indata.ntor+1): temp = temp + str(i-self.indata.ntor)+";"
			self.ui.TableArrays.setVerticalHeaderLabels(temp.split(";"))
			temp =''
			for i in range(self.indata.mpol): temp = temp + str(i)+";"
			self.ui.TableArrays.setHorizontalHeaderLabels(temp.split(";"))
			# Set Values
			val = getattr(self.indata,data_name)
			noffset = round((val.shape[1] - 1)/2)
			for i in range(2*self.indata.ntor+1):
				for m in range(self.indata.mpol):
					n1 = i - self.indata.ntor
					n2 = noffset + n1
					#i2 = noffset + i
					#i2 = noffset-self.indata.ntor + i -1
					self.ui.TableArrays.setItem(i,m, QTableWidgetItem(str(val[m,n2])))
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
				temp = getattr(self.indata,data_name+'_s')
				temp[col] = value
				setattr(self.indata,data_name+'_s',temp)
				#self.indata[data_name+'_s'][col]=value
				#set_module_var('vmec_input',data_name+'_s',self.indata[data_name+'_s'][:])
			else:
				temp = getattr(self.indata,data_name+'_f')
				temp[col] = value
				setattr(self.indata,data_name+'_f',temp)
				#self.indata[data_name+'_f'][col]=value
				#set_module_var('vmec_input',data_name+'_f',self.indata[data_name+'_f'][:])
		elif data_name in ['am','ac','ai']:
			temp = getattr(self.indata,data_name)
			temp[col] = value
			setattr(self.indata,data_name,temp)
			#self.indata[data_name][col]=value
			#set_module_var('vmec_input',data_name,self.indata[data_name][:])
		elif data_name == 'rbc' or data_name == 'zbs' or data_name == 'rbs' or data_name == 'zbc':
			temp = getattr(self.indata,data_name)
			row2 = 101 - self.indata.ntor + row
			temp[col][row2] = value
			setattr(self.indata,data_name,temp)
			#row2 = 101 - self.indata['ntor'] + row
			#self.indata[data_name][col][row2]=value
		self.DrawArrays()

	def DrawArrays(self):
		from libstell import libstell
		#from libstell import vmec
		# Determine type of array
		data_name = self.ui.ComboBoxArrays.currentText()
		data_name = data_name.lower()
		# Handle plots / values
		self.fig.clf()
		#self.ax = self.fig.add_axes([0, 0, 1, 1])
		self.ax = self.fig.add_subplot(111)
		self.canvas.show()
		self.vtkWidget.hide()
		s = np.linspace(0.0,1.0,99)
		#s = np.ndarray((99,1))
		f = np.zeros((99,1))
		#for i in range(99): s[i]=(i)/98.0
		if data_name[0:2] == 'am':
			for i,xx in enumerate(s):
				f[i] = self.indata.pmass(xx)/(pi*4E-7)
				self.ax.set_ylabel('Pressure')
		elif data_name[0:2] == 'ac':
			for i,xx in enumerate(s):
				f[i] = self.indata.pcurr(xx)
				self.ax.set_ylabel('Current')
		elif data_name[0:2] == 'ai':
			for i,xx in enumerate(s):
				f[i] = self.indata.piota(xx)
				self.ax.set_ylabel('Iota')
		elif data_name == 'rbc' or data_name == 'zbs' or data_name == 'rbs' or data_name == 'zbc':
			FOURIER_REP = libstell.FourierRep()
			mnmax = (2*self.indata.ntor+1)*(self.indata.mpol)
			nu = 4*self.indata.mpol
			nv = 4*self.indata.ntor
			if nu < 64: nu=64
			if nv < 32: nv=32
			mn = 0
			xn = np.ndarray((mnmax,1))
			xm = np.ndarray((mnmax,1))
			rmnc = np.ndarray((1,mnmax))
			zmns = np.ndarray((1,mnmax))
			rmns = np.ndarray((1,mnmax))
			zmnc = np.ndarray((1,mnmax))
			for i in range(2*self.indata.ntor+1):
				for j in range(self.indata.mpol):
					i2 = 101-self.indata.ntor + i
					xn[mn] = (i-self.indata.ntor)*self.indata.nfp
					xm[mn] = j
					rmnc[0,mn] = self.indata.rbc[j,i2]
					zmns[0,mn] = self.indata.zbs[j,i2]
					rmns[0,mn] = self.indata.rbs[j,i2]
					zmnc[0,mn] = self.indata.zbc[j,i2]
					mn = mn + 1
			theta = np.linspace([0],[np.pi*2],nu)
			zeta  = np.linspace([0],[np.pi/self.indata.nfp],nv)
			#theta = np.ndarray((nu,1))
			#zeta = np.ndarray((nv,1))
			#for j in range(nu): theta[j]=2*pi*j/(nu-1)
			#for j in range(nv): zeta[j]=pi*j/(nv*self.indata.nfp)
			r=FOURIER_REP.cfunct(theta,zeta,rmnc,xm,xn)
			z=FOURIER_REP.sfunct(theta,zeta,zmns,xm,xn)
			#self.ax = self.fig.add_subplot(111,projection='3d')
			self.canvas.hide()
			self.vtkWidget.show()
			self.plt_vmec.renderer.RemoveAllViewProps()
			self.ax = FOURIER_REP.isotoro(r,z,zeta,0,plot3D=self.plt_vmec,lclosev=False)
			self.canvas.draw()
			return
		fmax = max(f)
		if fmax == 0:
			fmax = 1
		#self.ax = self.fig.add_axes([0, 0, 1, 1])
		#self.ax = self.fig.add_subplot(111)
		self.ax.plot(s,f)
		self.ax.set_xlabel('Norm. Tor. Flux (s)')
		self.ax.set_aspect('auto')
		if data_name[3:6] == 'aux':
			x = getattr(self.indata,data_name+'_s')
			y = getattr(self.indata,data_name+'_f')
			self.ax.plot(x[x>=0],y[x>=0],'o')
		self.canvas.draw()

	def LoadOptimum(self):
		# Handles loading an indata file.
		w = QWidget()
		w.resize(320, 240)
		w.setWindowTitle("Hello World!")
		filename = QFileDialog.getOpenFileName(w, 'Open File', '.')
		w.destroy
		# Now read the file
		self.optimum.read_input(filename[0])
		self.indata.read_indata(filename[0])
		# Now update the UI
		if self.optimum.global_data.opt_type == 'one_iter':
			self.ui.ComboBoxOPTtype.setCurrentIndex(0)
		elif self.optimum.global_data.opt_type == 'lmdif':
			self.ui.ComboBoxOPTtype.setCurrentIndex(1)
		elif self.optimum.global_data.opt_type == 'lmdif_bounded':
			self.ui.ComboBoxOPTtype.setCurrentIndex(2)
		elif self.optimum.global_data.opt_type == 'gade':
			self.ui.ComboBoxOPTtype.setCurrentIndex(3)
		elif self.optimum.global_data.opt_type == 'pso':
			self.ui.ComboBoxOPTtype.setCurrentIndex(4)
		self.UpdateOPTtype()
		self.UpdateOPTVarsScalar()

	def UpdateOPTtype(self):
		self.ui.TableOPTtype.blockSignals(True)
		# Determine type of Optimization
		self.optimum.global_data.opt_type = self.ui.ComboBoxOPTtype.currentText()
		#self.optimum.global_data.opt_type = self.optimum.opt_type.upper()
		if self.optimum.global_data.opt_type == 'ONE_ITER':
			fields = ['noptimizers']
			self.ui.TableOPTtype.setRowCount(1)
		elif self.optimum.global_data.opt_type in ['LMDIF','LMDIF_BOUNDED']:
			self.ui.TableOPTtype.setRowCount(8)
			fields = ['nfunc_max','ftol','xtol','gtol','epsfcn','factor','mode','noptimizers']
		elif self.optimum.global_data.opt_type == 'GADE':
			self.ui.TableOPTtype.setRowCount(6)
			fields = ['nfunc_max','factor','cr_strategy','mode','npopulation','noptimizers']
		elif self.optimum.global_data.opt_type == 'PSO':
			self.ui.TableOPTtype.setRowCount(6)
			fields = ['nfunc_max','factor','cr_strategy','mode','npopulation','noptimizers']
		self.ui.TableOPTtype.setVerticalHeaderLabels(fields)
		for i,item in enumerate(fields):
			val = str(getattr(self.optimum.global_data,item))
		#	val = str(self.optimum[item])
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
		if field == 'opt_type':
			self.optimum.global_data.opt_type=item.text()
		elif field in ['noptimizers','nfunc_max','mode','cr_strategy','npopulation']:
			setattr(self.optimum.global_data,field,int(item.text()))
			#self.optimum[field]=int(item.text())
		else:
			setattr(self.optimum.global_data,field,float(item.text()))
			#self.optimum[field]=float(item.text())

	def StelVarsTab(self):
		# Handles updates to the tabs
		# Get page index
		item = self.ui.tabStelVars.currentIndex()
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
		for i,item in enumerate(['phiedge','pres_scale','curtor']):
			self.ui.TableOPTVarsScalar.setVerticalHeaderItem(i,QTableWidgetItem(item))
			vstate = getattr(self.indata,item)
			if item == 'pres_scale': item = 'pscale'
			lstate = getattr(self.optimum.var_data,'l'+item+'_opt')
			dstate = getattr(self.optimum.var_data,'d'+item+'_opt')
			minstate = getattr(self.optimum.var_data,item+'_min')
			maxstate = getattr(self.optimum.var_data,item+'_max')
			if lstate:
				self.ui.TableOPTVarsScalar.setItem(i,0,QTableWidgetItem('T'))
			else:
				self.ui.TableOPTVarsScalar.setItem(i,0,QTableWidgetItem('F'))
			self.ui.TableOPTVarsScalar.setItem(i,1,QTableWidgetItem(str(vstate)))
			self.ui.TableOPTVarsScalar.setItem(i,2,QTableWidgetItem(str(dstate)))
			self.ui.TableOPTVarsScalar.setItem(i,3,QTableWidgetItem(str(minstate)))
			self.ui.TableOPTVarsScalar.setItem(i,4,QTableWidgetItem(str(maxstate)))
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
		if field == 'pres_scale': field='pscale'
		# Column determins value
		field = field.lower()
		if col == 0:
			field = 'l'+field+'_opt'
			if item.text() == 'T':
				val = True
			else:
				val = False
		if col == 1:
			return
		elif col == 2:
			field = 'd'+field+'_opt'
			val = float(item.text())
		elif col == 3:
			field = field+'_min'
			val = float(item.text())
		elif col == 4:
			field = field+'_max'
			val = float(item.text())
		setattr(self.optimum.var_data,field,val)

	def UpdateOPTVarsProf(self):
		self.ui.TableOPTVarsProf.blockSignals(True)
		# Figure out which array to deal with
		data_name = self.ui.comboBoxStelVarsProfType.currentText()
		if data_name == 'Pressure (AM)':
			field = 'am'
		elif data_name == 'Current (AC)':
			field = 'ac'
		elif data_name == 'Iota (AI)':
			field = 'ai'
		if field in ['am','ac','ai']:
			nrows = 11;
		self.ui.TableOPTVarsProf.setRowCount(nrows)
		lstate = getattr(self.optimum.var_data,'l'+field+'_opt')
		vstate = getattr(self.indata,field)
		dstate = getattr(self.optimum.var_data,'d'+field+'_opt')
		minstate = getattr(self.optimum.var_data,field+'_min')
		maxstate = getattr(self.optimum.var_data,field+'_max')
		for i in range(nrows):
			self.ui.TableOPTVarsProf.setVerticalHeaderItem(i,QTableWidgetItem(str(field+'('+str(i)+')')))
			if lstate[i]:
				self.ui.TableOPTVarsProf.setItem(i,0,QTableWidgetItem('T'))
			else:
				self.ui.TableOPTVarsProf.setItem(i,0,QTableWidgetItem('F'))
			self.ui.TableOPTVarsProf.setItem(i,1,QTableWidgetItem(str(vstate[i])))
			self.ui.TableOPTVarsProf.setItem(i,2,QTableWidgetItem(str(dstate[i])))
			self.ui.TableOPTVarsProf.setItem(i,3,QTableWidgetItem(str(minstate[i])))
			self.ui.TableOPTVarsProf.setItem(i,4,QTableWidgetItem(str(maxstate[i])))
		self.ui.TableOPTVarsProf.blockSignals(False)

	def OPTVarsProf(self):
		# Handles changes in Profile table
		# Get the current profile
		data_name = self.ui.comboBoxStelVarsProfType.currentText()
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
		field = field.lower()
		if col == 0:
			field = 'l'+field+'_opt'
			if item.text() == 'T':
				val = True
			else:
				val = True
		if col == 1:
			return
		elif col == 2:
			field = 'd'+field+'_opt'
			val = float(item.text())
		elif col == 3:
			field = field+'_min'
			val = float(item.text())
		elif col == 4:
			field = field+'_max'
			val = float(item.text())
		temp = getattr(self.optimum.var_data,field)
		temp[row] = val
		setattr(self.optimum.var_data,field,temp)

	def UpdateOPTVarsExtcur(self):
		self.ui.TableOPTVarsExtcur.blockSignals(True)
		# See how many EXTCUR array vars there are
		if self.indata.lfreeb:
			nextcur = 0
			for i in range(99):
				if not(self.indata.extcur[i] == 0):
					nextcur = i
			self.ui.TableOPTVarsExtcur.setRowCount(nextcur)
		else:
			self.ui.TableOPTVarsExtcur.setRowCount(1)
			return
		field = 'extcur'
		lstate = getattr(self.optimum.var_data,'l'+field+'_opt')
		vstate = getattr(self.indata,field)
		dstate = getattr(self.optimum.var_data,'d'+field+'_opt')
		minstate = getattr(self.optimum.var_data,field+'_min')
		maxstate = getattr(self.optimum.var_data,field+'_max')
		for i in range(nextcur):
			if lstate[i]:
				self.ui.TableOPTVarsProf.setItem(i,0,QTableWidgetItem('T'))
			else:
				self.ui.TableOPTVarsProf.setItem(i,0,QTableWidgetItem('F'))
			self.ui.TableOPTVarsProf.setItem(i,1,QTableWidgetItem(str(vstate[i])))
			self.ui.TableOPTVarsProf.setItem(i,2,QTableWidgetItem(str(dstate[i])))
			self.ui.TableOPTVarsProf.setItem(i,3,QTableWidgetItem(str(minstate[i])))
			self.ui.TableOPTVarsProf.setItem(i,4,QTableWidgetItem(str(maxstate[i])))
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
		field = 'extcur'
		if col == 0:
			field = 'l'+field+'_opt'
			if item.text() == 'T':
				val = True
			else:
				val = False
		if col == 1:
			return
		elif col == 2:
			field = 'd'+field+'_opt'
			val = float(item.text())
		elif col == 3:
			field = field+'_min'
			val = float(item.text())
		elif col == 4:
			field = field+'_max'
			val = float(item.text())
		temp = getattr(self.optimum.var_data,field)
		temp[row] = val
		setattr(self.optimum.var_data,field,temp)

	def UpdateOPTVarsBound(self):
		self.ui.TableOPTVarsBound.blockSignals(True)
		# Figure out which array to deal with
		data_name = self.ui.comboBoxStelVarsBoundType.currentText()
		if data_name == 'VMEC':
			mmin = 0
			mmax = self.indata.mpol-1
			nmin = -self.indata.ntor
			nmax = self.indata.ntor
			field1 = 'bound'
			field2 = 'bound'
			field3 = 'bound'
			field4 = 'bound'
			msize = self.optimum.var_data.lbound_opt.shape[0]
			nsize = self.optimum.var_data.lbound_opt.shape[1]
			noffset = round((nsize-1)/2)
			moffset = 0
			self.ui.TableOPTVarsBound.setColumnCount(4)
			self.ui.TableOPTVarsBound.setHorizontalHeaderItem(0,QTableWidgetItem('Optimized'))
			self.ui.TableOPTVarsBound.setHorizontalHeaderItem(1,QTableWidgetItem('DValue'))
			self.ui.TableOPTVarsBound.setHorizontalHeaderItem(2,QTableWidgetItem('Minimum'))
			self.ui.TableOPTVarsBound.setHorizontalHeaderItem(3,QTableWidgetItem('Maximum'))
		elif data_name == 'Hirshman-Breslau':
			mmin = 0
			mmax = self.indata.mpol-1
			nmin = -self.indata.ntor
			nmax = self.indata.ntor
			field1 = 'rho'
			field2 = 'rho'
			field3 = 'bound'
			field4 = 'bound'
			msize = self.optimum.var_data.lrho_opt.shape[0]
			nsize = self.optimum.var_data.lrho_opt.shape[1]
			noffset = round((nsize-1)/2)
			moffset = 0
		elif data_name == 'Garabedian':
			mmin = -self.indata.mpol+1
			mmax = self.indata.mpol-1
			nmin = -self.indata.ntor
			nmax = self.indata.ntor
			field1 = 'deltamn'
			field2 = 'deltamn'
			field3 = 'delta'
			field4 = 'delta'
			msize = self.optimum.var_data.ldeltamn_opt.shape[0]
			nsize = self.optimum.var_data.ldeltamn_opt.shape[1]
			noffset = round((nsize-1)/2)
			moffset = round((msize-1)/2)
		# Get var data
		lstate = getattr(self.optimum.var_data,'l'+field1+'_opt')
		dstate = getattr(self.optimum.var_data,'d'+field2+'_opt')
		minstate = getattr(self.optimum.var_data,field3+'_min')
		maxstate = getattr(self.optimum.var_data,field4+'_max')
		xm    = []
		xn    = []
		ltemp = []
		dtemp = []
		mintemp = []
		maxtemp = []
		# Create list of values with (n,m) as row headers
		for j in range(msize):
			for i in range(nsize):
				ng = i - noffset
				mg = j - moffset
				if ng >= nmin and ng <= nmax and \
					mg >= mmin and mg <= mmax:
					xm.append(mg)
					xn.append(ng)
					ltemp.append(lstate[j,i])
					dtemp.append(dstate[j,i])
					mintemp.append(minstate[j,i])
					maxtemp.append(maxstate[j,i])
		# Setup the Array
		j = len(xm)
		self.ui.TableOPTVarsBound.setRowCount(j)
		for i in range(j):
			head_string = '('+str(xn[i])+','+str(xm[i])+')'
			self.ui.TableOPTVarsBound.setVerticalHeaderItem(i,QTableWidgetItem(head_string))
			if ltemp[i]:
				self.ui.TableOPTVarsBound.setItem(i,0,QTableWidgetItem('T'))
			else:
				self.ui.TableOPTVarsBound.setItem(i,0,QTableWidgetItem('F'))
			self.ui.TableOPTVarsBound.setItem(i,1,QTableWidgetItem(str(dtemp[i])))
			self.ui.TableOPTVarsBound.setItem(i,2,QTableWidgetItem(str(mintemp[i])))
			self.ui.TableOPTVarsBound.setItem(i,3,QTableWidgetItem(str(maxtemp[i])))
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
			field1 = 'bound'
			field2 = 'bound'
			field3 = 'bound'
			field4 = 'bound'
			dex1 = n+100
			dex2 = m
		# Column determins value
		if col == 0:
			field = 'l'+field1+'_opt'
			if item.text() == 'T':
				val = True
			else:
				val = False
		if col == 1:
			return
		elif col == 2:
			field = 'd'+field2+'_opt'
			val = float(item.text())
		elif col == 3:
			field = field3+'_min'
			val = float(item.text())
		elif col == 4:
			field = field4+'_max'
			val = float(item.text())
		temp = getattr(self.optimum.var_data,field)
		temp[dex1,dex2] = val
		setattr(self.optimum.var_data,field,temp)

	def LoadSTELLOPT(self):
		# Handles loading an stellopt file.
		w = QWidget()
		w.resize(320, 240)
		w.setWindowTitle("Load STELLOPT Output")
		filename = QFileDialog.getOpenFileName(w, 'Open File', '.','STELLOPT (stellopt.*)')
		w.destroy
		# Read the file
		self.stel_data.read_stellopt_output(filename[0])
		self.optplot_list = ['ASPECT','BETA','CURTOR','EXTCUR','SEPARATRIX',\
					'PHIEDGE','RBTOR','R0','Z0','VOLUME','WP','KAPPA',\
					'B_PROBES','FARADAY','FLUXLOOPS','SEGROG','MSE',\
					'NE','NELINE','TE','TELINE','TI','TILINE','ZEFFLINE',\
					'XICS','XICS_BRIGHT','XICS_W3','XICS_V','SXR','VPHI','VACIOTA',\
					'IOTA','BALLOON','BOOTSTRAP','DKES','DKES_ERDIFF',\
					'HELICITY','HELICITY_FULL',\
					'KINK','ORBIT','JDOTB','J_STAR','NEO','TXPORT','ECEREFLECT',\
					'S11','S12','S21','S22','MAGWELL',\
					'CURVATURE_KERT','CURVATURE_P2','QUASIISO']
		self.ui.ComboBoxOPTplot_type.clear()
		self.ui.ComboBoxOPTplot_type.addItem('Chi-Squared')
		# Handle Chisquared plots
		for name in self.optplot_list:
			for item in vars(self.stel_data).keys():
				if (name+'_TARGET' == item):
					self.ui.ComboBoxOPTplot_type.addItem(name)
		# Handle Special Plots
		self.ui.ComboBoxOPTplot_type.addItem('-----SPECIAL-----')
		for name in ['BALLOON','KINK','ORBIT','NEO','HELICITY','HELICITY_FULL',\
					'TXPORT','B_PROBES','FLUXLOOPS','SEGROG',\
					'NELINE','TELINE','TILINE','ZEFFLINE',\
					'XICS','XICS_BRIGHT','XICS_W3','XICS_V',\
					'S11','S12','S21','S22','MAGWELL','VACIOTA',\
					'CURVATURE_KERT','CURVATURE_P2',\
					'ECEREFLECT','SXR','IOTA','PRESS','PRESSPRIME'\
					'VISBREMLINE','QUASIISO']:
			for item in vars(self.stel_data).keys():
				if (name+'_TARGET' == item):
					self.ui.ComboBoxOPTplot_type.addItem(name+'_evolution')
		for name in ['NE','TE','TI','MSE']:
			for item in vars(self.stel_data).keys():
				if (name+'_TARGET' == item):
					self.ui.ComboBoxOPTplot_type.addItem(name+'_evolution')
					self.ui.ComboBoxOPTplot_type.addItem(name+'_evolution_R')
					self.ui.ComboBoxOPTplot_type.addItem(name+'_evolution_Z')
		if 'DKES_TARGET' in vars(self.stel_data).keys():
			self.ui.ComboBoxOPTplot_type.addItem('DKES_L11')
			self.ui.ComboBoxOPTplot_type.addItem('DKES_L31')
			self.ui.ComboBoxOPTplot_type.addItem('DKES_L33')
		# Handle Wout Comparrison Plots
		self.workdir,ext = filename[0].split('stellopt.',1)
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
			self.ui.ComboBoxOPTplot_type.addItem('Mercier')
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
		niter = len(self.stel_data.ITER)
		#self.fig.delaxes(self.ax)
		#self.ax2 = self.fig2.add_subplot(111)
		self.ax2 = self.fig2.add_axes([0.2,0.2,0.7,0.7])
		if (plot_name == 'Chi-Squared'):
			chisq = ((self.stel_data.TARGETS - self.stel_data.VALS)/self.stel_data.SIGMAS)**2
			self.ax2.semilogy(self.stel_data.ITER,np.sum(chisq,axis=1),'ok',label='Chisq Total')
			self.ax2.set_xlabel('Iteration')
			self.ax2.set_ylabel('Chi-Squared')
			self.ax2.set_title('Chi-Sqaured')
			#self.ax2.set_yscale('log')
			for name in self.optplot_list:
				if name+'_CHISQ' in vars(self.stel_data).keys():
					chisq_temp = getattr(self.stel_data,name+'_CHISQ')
					#chisq_temp = self.stel_data[name+'_chisq']
					n = chisq_temp.shape;
					if (len(chisq_temp.shape) == 1):
						if n[0] > len(self.stel_data.ITER):
							chisq_temp = np.sum(chisq_temp,axis=0)
					elif len(chisq_temp.shape) > 1:
						chisq_temp = np.sum(chisq_temp,axis=1)
					self.ax2.plot(self.stel_data.ITER,chisq_temp,'o',fillstyle='none',label=name)
			self.ax2.legend()
		elif (plot_name in self.optplot_list):
			f = getattr(self.stel_data,plot_name+'_CHISQ')
			#f = self.stel_data[plot_name+'_chisq']
			n = f.shape
			if (len(n)==0):
				# Single Time slice Single point
				n=0
			elif (len(n)==1):
				# Could be either mutli-time or single time
				if len(self.stel_data.ITER) == n[0]:
					# Mutl-time single point
					n=0
				else:
					# Multi-channel single time
					f = np.sum(f,axis=0)
			else:
				# Multiple Time slices
				f = np.sum(f,axis=1)
			self.ax2.semilogy(self.stel_data.ITER,f,'ok',fillstyle='none')
			self.ax2.set_xlabel('Iteration')
			self.ax2.set_ylabel('Chi-Squared')
			self.ax2.set_title(plot_name+' Chi-Sqaured')
			#self.ax2.set_yscale('log',basey=10)
		elif (plot_name == 'BALLOON_evolution'):
			x = self.stel_data.BALLOON_K
			y = self.stel_data.BALLOON_BALLOON_GRATE
			t = self.stel_data.BALLOON_TARGET
			d = self.stel_data.BALLOON_SIGMA
			self.ax2.errorbar(x[0,:],t[0,:],yerr=d[0,:],fmt='ok',fillstyle='none',label='Target')
			self.ax2.plot(x[0,:],y[0,:],'o',fillstyle='none',label='Initial',color='red')
			for i in range(1,niter-1,1):
				self.ax2.plot(x[i,:],y[i,:],'.k',fillstyle='none')
			self.ax2.plot(x[niter-1,:],y[niter-1,:],'o',fillstyle='none',label='Final',color='green')
			self.ax2.set_xlabel('Radial Grid')
			self.ax2.set_ylabel('Growth Rate')
			self.ax2.set_title('COBRA Ballooning Stability (<0 Stable)')
			self.ax2.legend()
		elif (plot_name == 'TXPORT_evolution'):
			x = self.stel_data.TXPORT_S
			y = self.stel_data.TXPORT_VAL
			t = self.stel_data.TXPORT_TARGET
			d = self.stel_data.TXPORT_SIGMA
			self.ax2.errorbar(x[0,:],t[0,:],yerr=d[0,:],fmt='ok',fillstyle='none',label='Target')
			self.ax2.plot(x[0,:],y[0,:],'o',fillstyle='none',label='Initial',color='red')
			for i in range(1,niter-1,1):
				self.ax2.plot(x[i,:],y[i,:],'.k',fillstyle='none')
			self.ax2.plot(x[niter-1,:],y[niter-1,:],'o',fillstyle='none',label='Final',color='green')
			self.ax2.set_xlabel('Normalized Flux')
			self.ax2.set_ylabel('Proxy Function')
			self.ax2.set_title('Turbulent Transport Proxy')
			self.ax2.set_xlim((0,1))
		elif (plot_name == 'ORBIT_evolution'):
			x = self.stel_data.ORBIT_S
			y = self.stel_data.ORBIT_EQUIL
			t = self.stel_data.ORBIT_TARGET
			d = self.stel_data.ORBIT_SIGMA
			self.ax2.errorbar(x[0,:],t[0,:],yerr=d[0,:],fmt='ok',fillstyle='none',label='Target')
			self.ax2.plot(x[0,:],y[0,:],'o',fillstyle='none',label='Initial',color='red')
			for i in range(1,niter-1,1):
				self.ax2.plot(x[i,:],y[i,:],'.k',fillstyle='none')
			self.ax2.plot(x[niter-1,:],y[niter-1,:],'o',fillstyle='none',label='Final',color='green')
			self.ax2.set_xlabel('Normalized Flux')
			self.ax2.set_ylabel('Orbit Losses')
			self.ax2.set_title('Gyro Particle Losses')
			self.ax2.set_xlim((0,1))
		elif (plot_name == 'NEO_evolution'):
			x = self.stel_data.NEO_K
			y = self.stel_data.NEO_EPS_EFF32
			t = self.stel_data.NEO_TARGET
			d = self.stel_data.NEO_SIGMA
			self.ax2.errorbar(x[0,:],t[0,:],yerr=d[0,:],fmt='ok',fillstyle='none',label='Target')
			self.ax2.plot(x[0,:],y[0,:],'o',fillstyle='none',label='Initial',color='red')
			for i in range(1,niter-1,1):
				self.ax2.plot(x[i,:],y[i,:],'.k',fillstyle='none')
			self.ax2.plot(x[niter-1,:],y[niter-1,:],'o',fillstyle='none',label='Final',color='green')
			self.ax2.set_xlabel('Radial Grid')
			self.ax2.set_ylabel('Epsilon Effective')
			self.ax2.set_title('Neoclassical Helical Ripple (NEO)')
		elif ('DKES_L' in plot_name):
			# Get L type
			if plot_name == 'DKES_L11':
				Lm = self.stel_data.DKES_L11m
				Lp = self.stel_data.DKES_L11p
				txt_type = 'L11'
			elif plot_name == 'DKES_L31':
				Lm = self.stel_data.DKES_L31m
				Lp = self.stel_data.DKES_L31p
				txt_type = 'L31'
			elif plot_name == 'DKES_L33':
				Lm = self.stel_data.DKES_L33m
				Lp = self.stel_data.DKES_L33p
				txt_type = 'L33'
			# We need to sort stuff out
			s  = self.stel_data.DKES_S
			er = self.stel_data.DKES_ER
			nu = self.stel_data.DKES_NU
			s_list  = np.unique(s)
			er_list = np.unique(er)
			nu_list = np.unique(nu)
			ns  = len(s_list)
			ner = len(er_list)
			nnu = len(nu_list)
			s3d  = np.zeros((ns,ner,nnu))
			er3d = np.zeros((ns,ner,nnu))
			nu3d = np.zeros((ns,ner,nnu))
			Lval  = np.zeros((ns,ner,nnu,2))
			for i in range(ns):
				for j in range(ner):
					for k in range(nnu):
						sdex = s == s_list[i]
						edex = er == er_list[j]
						ndex = nu == nu_list[k]
						dex = np.logical_and(sdex,edex)
						dex = np.logical_and(dex,ndex)
						s3d[i,j,k] = np.squeeze(s[dex])
						er3d[i,j,k] = np.squeeze(er[dex])
						nu3d[i,j,k] = np.squeeze(nu[dex])
						Lval[i,j,k,0] = np.squeeze(Lm[dex])
						Lval[i,j,k,1] = np.squeeze(Lp[dex])
			for i in range(ns):
				for j in range(ner):
					x = np.squeeze(nu3d[i,j,:])
					ym = np.squeeze(Lval[i,j,:,0])
					yp = np.squeeze(Lval[i,j,:,1])
					self.ax2.fill_between(x,ym,yp,alpha=0.2)
			self.ax2.set_xlabel('Collisionality nu*')
			self.ax2.set_ylabel(txt_type)
			self.ax2.set_title("DKES Coefficient "+txt_type)
			self.ax2.set_yscale('log')
		elif (plot_name == 'HELICITY_FULL_evolution'):
			x = self.stel_data.HELICITY_FULL_K
			y = self.stel_data.HELICITY_FULL_VAL
			t = self.stel_data.HELICITY_FULL_TARGET
			d = self.stel_data.HELICITY_FULL_SIGMA
			self.ax2.errorbar(x[0,:],t[0,:],yerr=d[0,:],fmt='ok',fillstyle='none',label='Target')
			self.ax2.plot(x[0,:],y[0,:],'o',fillstyle='none',label='Initial',color='red')
			for i in range(1,niter-1,1):
				self.ax2.plot(x[i,:],y[i,:],'.k',fillstyle='none')
			self.ax2.plot(x[niter-1,:],y[niter-1,:],'o',fillstyle='none',label='Final',color='green')
			self.ax2.set_ylabel('Helicity')
			self.ax2.set_title('Boozer Spectrum Helicity')
		elif (plot_name == 'MAGWELL_evolution'):
			x = self.stel_data.MAGWELL_k
			y = self.stel_data.MAGWELL_MAGWELL
			t = self.stel_data.MAGWELL_TARGET
			d = self.stel_data.MAGWELL_SIGMA
			self.ax2.errorbar(x[0,:],t[0,:],yerr=d[0,:],fmt='ok',fillstyle='none',label='Target')
			self.ax2.plot(x[0,:],y[0,:],'o',fillstyle='none',label='Initial',color='red')
			for i in range(1,niter-1,1):
				self.ax2.plot(x[i,:],y[i,:],'.k',fillstyle='none')
			self.ax2.plot(x[niter-1,:],y[niter-1,:],'o',fillstyle='none',label='Final',color='green')
			self.ax2.set_xlabel('Radial Grid')
			self.ax2.set_ylabel('Magnetic Well')
			self.ax2.set_title('Magnetic Well Evolution  (>0 Well)')
		elif (plot_name == 'QUASIISO_evolution'):
			x = self.stel_data.QUASIISO_K
			y = self.stel_data.QUASIISO_VAL
			t = self.stel_data.QUASIISO_TARGET
			d = self.stel_data.QUASIISO_SIGMA
			self.ax2.errorbar(x[0,:],t[0,:],yerr=d[0,:],fmt='ok',fillstyle='none',label='Target')
			self.ax2.plot(x[0,:],y[0,:],'o',fillstyle='none',label='Initial',color='red')
			for i in range(1,niter-1,1):
				self.ax2.plot(x[i,:],y[i,:],'.k',fillstyle='none')
			self.ax2.plot(x[niter-1,:],y[niter-1,:],'o',fillstyle='none',label='Final',color='green')
			self.ax2.set_xlabel('Radial Grid')
			self.ax2.set_ylabel('QI Metric')
			self.ax2.set_title('Quasi-isodynamic Metric')
		elif (plot_name == 'CURVATURE_P2'):
			#x = self.stel_data.TXPORT_S
			#y = self.stel_data.CURVATURE_P2_P2
			#self.ax2.plot(x[0,:],y[0,:],'o',fillstyle='none',label='Initial',color='red')
			#for i in range(1,niter-1,1):
			#	self.ax2.plot(x[i,:],y[i,:],'.k',fillstyle='none')
			#self.ax2.plot(x[niter-1,:],self.y[niter-1,:],'o',fillstyle='none',label='Final',color='green')
			self.ax2.plot(self.stel_data.CURVATURE_P2_P2.T,'o',fillstyle='none')
			self.ax2.set_xlabel('Radial Grid')
			self.ax2.set_ylabel('Magnetic Well')
			self.ax2.set_title('Magnetic Well Evolution')
		elif (plot_name == 'B_PROBE_evolution'):
			n=self.stel_data['B_PROBE_target'].shape
			x = np.ndarray((niter,n[1]))
			y = self.stel_data.BPROBES_equil
			t = self.stel_data.BPROBES_TARGET
			d = self.stel_data.BPROBES_SIGMA
			self.ax2.errorbar(x[0,:],t[0,:],yerr=d[0,:],fmt='ok',fillstyle='none',label='Target')
			self.ax2.plot(x[0,:],y[0,:],'o',fillstyle='none',label='Initial',color='red')
			for i in range(1,niter-1,1):
				self.ax2.plot(x[i,:],y[i,:],'.k',fillstyle='none')
			self.ax2.plot(x[niter-1,:],y[niter-1,:],'o',fillstyle='none',label='Final',color='green')
			self.ax2.set_xlabel('B-Probe')
			self.ax2.set_ylabel('Signal')
			self.ax2.set_title('B-Probe Reconstruction')
		elif (plot_name == 'FLUXLOOPS_evolution'):
			n=self.stel_data.FLUXLOOPS_TARGET.shape
			y = self.stel_data.FLUXLOOPS_TARGET.T
			s = self.stel_data.FLUXLOOPS_SIGMA.T
			e = self.stel_data.FLUXLOOPS_EQUIL.T
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
			n=self.stel_data.SEGROG_TARGET.shape
			y = self.stel_data.SEGROG_TARGET.T
			s = self.stel_data.SEGROG_SIGMA.T
			e = self.stel_data.SEGROG_EQUIL.T
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
			n=self.stel_data.ECEREFLECT_TARGET.shape
			y = self.stel_data.ECEREFLECT_TARGET.T
			s = self.stel_data.ECEREFLECT_SIGMA.T
			e = self.stel_data.ECEREFLECT_EQUIL.T
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
			self.ax2.plot(self.stel_data.KINK_EQUIL,fmt='ok',fillstyle='none')
			self.ax2.plot(self.stel_data.KINK_TARGET,fmt='k')
			self.ax2.plot(self.stel_data.KINK_TARGET+self.stel_data.KINK_SIGMA,fmt='k')
			self.ax2.plot(self.stel_data.KINK_TARGET-self.stel_data.KINK_SIGMA,fmt='k')
			self.ax2.set_xlabel('Iteration')
			self.ax2.set_ylabel('???')
			self.ax2.set_title('?????KINK Evolution????')
		elif (plot_name == 'NE_evolution'):
			x = self.stel_data.NE_S.T
			y = self.stel_data.NE_TARGET.T
			s = self.stel_data.NE_SIGMA.T
			e = self.stel_data.NE_VAL.T
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
			x = self.stel_data.NE_R.T
			y = self.stel_data.NE_TARGET.T
			s = self.stel_data.NE_SIGMA.T
			e = self.stel_data.NE_VAL.T
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
			x = self.stel_data.NE_Z.T
			y = self.stel_data.NE_TARGET.T
			s = self.stel_data.NE_SIGMA.T
			e = self.stel_data.NE_EQUIL.T
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
			x = self.stel_data.TE_S.T
			y = self.stel_data.TE_TARGET.T
			s = self.stel_data.TE_SIGMA.T
			e = self.stel_data.TE_VAL.T
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
			x = self.stel_data.TE_R.T
			y = self.stel_data.TE_TARGET.T
			s = self.stel_data.TE_SIGMA.T
			e = self.stel_data.TE_VAL.T
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
			x = self.stel_data.TE_Z.T
			y = self.stel_data.TE_TARGET.T
			s = self.stel_data.TE_SIGMA.T
			e = self.stel_data.TE_VAL.T
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
			x = self.stel_data.TI_S.T
			y = self.stel_data.TI_TARGET.T
			s = self.stel_data.TI_SIGMA.T
			e = self.stel_data.TI_VAL.T
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
			x = self.stel_data.TI_R.T
			y = self.stel_data.TI_TARGET.T
			s = self.stel_data.TI_SIGMA.T
			e = self.stel_data.TI_VAL.T
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
			x = self.stel_data.TI_Z.T
			y = self.stel_data.TI_TARGET.T
			s = self.stel_data.TI_SIGMA.T
			e = self.stel_data.TI_VAL.T
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
			x = self.stel_data.MSE_S.T
			y = self.stel_data.MSE_TARGET.T
			s = self.stel_data.MSE_SIGMA.T
			e = self.stel_data.MSE_VAL.T
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
			x = self.stel_data.MSE_R.T
			y = self.stel_data.MSE_TARGET.T
			s = self.stel_data.MSE_SIGMA.T
			e = self.stel_data.MSE_VAL.T
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
			x = self.stel_data.MSE_Z.T
			y = self.stel_data.MSE_TARGET.T
			s = self.stel_data.MSE_SIGMA.T
			e = self.stel_data.MSE_VAL.T
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
			x = self.stel_data.IOTA_S.T
			y = self.stel_data.IOTA_TARGET.T
			s = self.stel_data.IOTA_SIGMA.T
			e = self.stel_data.IOTA_VAL.T
			n = y.shape
			if len(x.shape)>1:
				dl = n[1]
				self.ax2.errorbar(x[:,0],y[:,0],s[:,0],fmt='sk',fillstyle='none')
				for l in range(dl):
					self.ax2.plot(x[:,l],e[:,l],'o',fillstyle='none',color=_plt.cm.brg(l/max(dl-1,1.0)))
			else:
				self.ax2.errorbar(x[:],y[:],s[:],fmt='sk',fillstyle='none')
				self.ax2.plot(x[:],e[:],'o',fillstyle='none',color='g')
			self.ax2.set_xlabel('Normalized Flux')
			self.ax2.set_ylabel('Iota')
			self.ax2.set_title('Rotational Transform')
		elif (plot_name == 'NELINE_evolution'):
			y=self.stel_data.NELINE_TARGET.T
			s=self.stel_data.NELINE_SIGMA.T
			e = self.stel_data.NELINE_VAL.T
			n = y.shape
			if (len(n)==0):
				# Single Time slice Single point
				x=np.ndarray((1,1))*0+1
				self.ax2.errorbar(x,y,s,fmt='sk',fillstyle='none')
				self.ax2.plot(x,e,'og',fillstyle='none')
			elif (len(n)==1):
				# Could be either mutli-time or single time
				if len(self.stel_data.ITER) == n[0]:
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
			y=self.stel_data.TELINE_TARGET.T
			s=self.stel_data.TELINE_SIGMA.T
			e = self.stel_data.TELINE_VAL.T
			n = y.shape
			dl = n[0]
			if (len(n)==0):
				# Single Time slice Single point
				x=np.ndarray((1,1))*0+1
				self.ax2.errorbar(x,y,s,fmt='sk',fillstyle='none')
				self.ax2.plot(x,e,'og',fillstyle='none')
			elif (len(n)==1):
				# Could be either mutli-time or single time
				if len(self.stel_data.ITER) == n[0]:
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
			y=self.stel_data.TILINE_TARGET.T
			s=self.stel_data.TILINE_SIGMA.T
			e = self.stel_data.TILINE_VAL.T
			n = y.shape
			if (len(n)==0):
				# Single Time slice Single point
				x=np.ndarray((1,1))*0+1
				self.ax2.errorbar(x,y,s,fmt='sk',fillstyle='none')
				self.ax2.plot(x,e,'og',fillstyle='none')
			elif (len(n)==1):
				# Could be either mutli-time or single time
				if len(self.stel_data.ITER) == n[0]:
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
		elif (plot_name == 'ZEFFLINE_evolution'):
			y=self.stel_data.ZEFFLINE_TARGET.T
			s=self.stel_data.ZEFFLINE_SIGMA.T
			e = self.stel_data.ZEFFLINE_VAL.T
			n = y.shape
			if (len(n)==0):
				# Single Time slice Single point
				x=np.ndarray((1,1))*0+1
				self.ax2.errorbar(x,y,s,fmt='sk',fillstyle='none')
				self.ax2.plot(x,e,'og',fillstyle='none')
			elif (len(n)==1):
				# Could be either mutli-time or single time
				if len(self.stel_data.ITER) == n[0]:
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
			self.ax2.set_title('Line-Int. Z Effective')
		elif (plot_name == 'XICS_evolution'):
			y=self.stel_data.XICS_TARGET.T
			s=self.stel_data.XICS_SIGMA.T
			e = self.stel_data.XICS_VAL.T
			n = y.shape
			if (len(n)==0):
				# Single Time slice Single point
				x=np.ndarray((1,1))*0+1
				self.ax2.errorbar(x,y,s,fmt='sk',fillstyle='none')
				self.ax2.plot(x,e,'og',fillstyle='none')
			elif (len(n)==1):
				# Could be either mutli-time or single time
				if len(self.stel_data.ITER) == n[0]:
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
			y=self.stel_data.XICS_BRIGHT_TARGET.T
			s=self.stel_data.XICS_BRIGHT_SIGMA.T
			e = self.stel_data.XICS_BRIGHT_VAL.T
			n = y.shape
			if (len(n)==0):
				# Single Time slice Single point
				x=np.ndarray((1,1))*0+1
				self.ax2.errorbar(x,y,s,fmt='sk',fillstyle='none')
				self.ax2.plot(x,e,'og',fillstyle='none')
			elif (len(n)==1):
				# Could be either mutli-time or single time
				if len(self.stel_data.ITER) == n[0]:
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
			y=self.stel_data.XICS_W3_TARGET.T
			s=self.stel_data.XICS_W3_SIGMA.T
			e = self.stel_data.XICS_W3_VAL.T
			n = y.shape
			if (len(n)==0):
				# Single Time slice Single point
				x=np.ndarray((1,1))*0+1
				self.ax2.errorbar(x,y,s,fmt='sk',fillstyle='none')
				self.ax2.plot(x,e,'og',fillstyle='none')
			elif (len(n)==1):
				# Could be either mutli-time or single time
				if len(self.stel_data.ITER) == n[0]:
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
			y=self.stel_data.XICS_V_TARGET.T
			s=self.stel_data.XICS_V_SIGMA.T
			e = self.stel_data.XICS_V_VAL.T
			n = y.shape
			if (len(n)==0):
				# Single Time slice Single point
				x=np.ndarray((1,1))*0+1
				self.ax2.errorbar(x,y,s,fmt='sk',fillstyle='none')
				self.ax2.plot(x,e,'og',fillstyle='none')
			elif (len(n)==1):
				# Could be either mutli-time or single time
				if len(self.stel_data.ITER) == n[0]:
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
			vmec_data = vmec.VMEC()
			l=0
			dl = len(self.wout_files)-1
			if dl == 0 : dl = 1 
			for string in self.wout_files:
				if 'wout' in string:
					vmec_data.read_wout(self.workdir+string)
					ns = vmec_data.ns
					nflux = np.ndarray((ns,1))
					for j in range(ns): nflux[j]=j/(ns-1)
					self.ax2.plot(nflux,vmec_data.presf/1000,color=_plt.cm.brg(l/dl))
					l=l+1
			self.ax2.set_xlabel('Norm Tor. Flux (s)')
			self.ax2.set_ylabel('Pressure [kPa]')
			self.ax2.set_title('VMEC Pressure Evolution')
			self.ax2.set_xlim((0,1))
		elif (plot_name == 'I-prime'):
			vmec_data = vmec.VMEC()
			l=0
			dl = len(self.wout_files)-1
			if dl == 0 : dl = 1 
			for string in self.wout_files:
				if 'wout' in string:
					vmec_data.read_wout(self.workdir+string)
					ns = vmec_data.ns
					nflux = np.ndarray((ns,1))
					for j in range(ns): nflux[j]=j/(ns-1)
					self.ax2.plot(nflux,vmec_data.jcurv/1000,color=_plt.cm.brg(l/dl))
					l=l+1
			self.ax2.set_xlabel('Norm Tor. Flux (s)')
			self.ax2.set_ylabel('Current Density dI/ds [kA]')
			self.ax2.set_title('VMEC Current Density Evolution')
			self.ax2.set_xlim((0,1))
		elif (plot_name == 'Iota'):
			vmec_data = vmec.VMEC()
			l=0
			dl = len(self.wout_files)-1
			for string in self.wout_files:
				if 'wout' in string:
					vmec_data.read_wout(self.workdir+string)
					ns = vmec_data.ns
					nflux = np.ndarray((ns,1))
					for j in range(ns): nflux[j]=j/(ns-1)
					self.ax2.plot(nflux,vmec_data.iotaf,color=_plt.cm.brg(l/dl))
					l=l+1
			self.ax2.set_xlabel('Norm Tor. Flux (s)')
			self.ax2.set_ylabel(r'Rotational Transform $\iota$')
			self.ax2.set_title('VMEC Rotational Transform Evolution')
			self.ax2.set_xlim((0,1))
		elif (plot_name == 'q-prof'):
			vmec_data = vmec.VMEC()
			l=0
			dl = len(self.wout_files)
			for string in self.wout_files:
				if 'wout' in string:
					vmec_data.read_wout(self.workdir+string)
					ns = vmec_data.ns
					nflux = np.ndarray((ns,1))
					for j in range(ns): nflux[j]=j/(ns-1)
					self.ax2.plot(nflux,1.0/vmec_data.iotaf,color=_plt.cm.brg(l/dl))
					l=l+1
			self.ax2.set_xlabel('Norm Tor. Flux (s)')
			self.ax2.set_ylabel('Safety Factor q')
			self.ax2.set_title('VMEC Safety Factor Evolution')
			self.ax2.set_xlim((0,1))
		elif (plot_name == 'Current'):
			vmec_data = vmec.VMEC()
			l=0
			dl = len(self.wout_files)-1
			for string in self.wout_files:
				if 'wout' in string:
					vmec_data.read_wout(self.workdir+string)
					ns = vmec_data.ns
					nflux = np.ndarray((ns,1))
					for j in range(ns): nflux[j]=j/(ns-1)
					self.ax2.plot(nflux,np.cumsum(vmec_data.jcurv),color=_plt.cm.brg(l/dl))
					l=l+1
			self.ax2.set_xlabel('Norm Tor. Flux (s)')
			self.ax2.set_ylabel('Current [kA]')
			self.ax2.set_title('VMEC Current Evolution')
			self.ax2.set_xlim((0,1))
		elif (plot_name == '<j*B>'):
			vmec_data = vmec.VMEC()
			l=0
			dl = len(self.wout_files)-1
			for string in self.wout_files:
				if 'wout' in string:
					vmec_data.read_wout(self.workdir+string)
					ns = vmec_data.ns
					nflux = np.ndarray((ns,1))
					for j in range(ns): nflux[j]=j/(ns-1)
					self.ax2.plot(nflux,vmec_data.jdotb,color=_plt.cm.brg(l/dl))
					l=l+1
			self.ax2.set_xlabel('Norm Tor. Flux (s)')
			self.ax2.set_ylabel('<j.B>')
			self.ax2.set_title('VMEC <j.B> Evolution')
			self.ax2.set_xlim((0,1))
		elif (plot_name == 'Mercier'):
			vmec_data = vmec.VMEC()
			l=0
			dl = len(self.wout_files)-1
			for string in self.wout_files:
				if 'wout' in string:
					vmec_data.read_wout(self.workdir+string)
					ns = vmec_data.ns
					nflux = np.ndarray((ns,1))
					for j in range(ns): nflux[j]=j/(ns-1)
					self.ax2.plot(nflux,vmec_data.dmerc,color=_plt.cm.brg(l/dl))
					l=l+1
			self.ax2.set_xlabel('Norm Tor. Flux (s)')
			self.ax2.set_ylabel('[Arb]')
			self.ax2.set_title('Mercier Stability (>0 Stable)')
			self.ax2.set_xlim((0,1))
		elif (plot_name == 'Flux0'):
			vmec_data = vmec.VMEC()
			l=0
			dl = len(self.wout_files)-1
			for string in self.wout_files:
				if 'wout' in string:
					vmec_data.read_wout(self.workdir+string)
					ns = vmec_data.ns
					nu = 256
					nv = 1
					nflux = np.ndarray((ns,1))
					theta = np.ndarray((nu,1))
					zeta = np.ndarray((nv,1))
					for j in range(ns): nflux[j]=j/(ns-1)
					for j in range(nu): theta[j]=2*pi*j/(nu-1)
					zeta[0]=0;
					r=vmec_data.cfunct(theta,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn)
					z=vmec_data.sfunct(theta,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn)
					self.ax2.plot(r[ns-1,:,0],z[ns-1,:,0],color=_plt.cm.brg(l/dl))
					self.ax2.plot(r[0,0,0],z[0,0,0],'+',color=_plt.cm.brg(l/dl))
					l=l+1
			self.ax2.set_xlabel('R [m]')
			self.ax2.set_ylabel('Z [m]')
			self.ax2.set_title('VMEC Flux Surface Evolution (phi=0)')
			self.ax2.set_aspect('equal')
		elif (plot_name == 'FluxPI'):
			vmec_data = vmec.VMEC()
			l=0
			dl = len(self.wout_files)-1
			if dl == 0 : dl = 1 
			for string in self.wout_files:
				if 'wout' in string:
					vmec_data.read_wout(self.workdir+string)
					ns = vmec_data.ns
					nu = 256
					nv = 1
					nflux = np.ndarray((ns,1))
					theta = np.ndarray((nu,1))
					zeta = np.ndarray((nv,1))
					for j in range(ns): nflux[j]=j/(ns-1)
					for j in range(nu): theta[j]=2*pi*j/(nu-1)
					zeta[0]=pi/vmec_data.nfp;
					r=vmec_data.cfunct(theta,zeta,vmec_data.rmnc,vmec_data.xm,vmec_data.xn)
					z=vmec_data.sfunct(theta,zeta,vmec_data.zmns,vmec_data.xm,vmec_data.xn)
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

	def PlotSTELLOPT(self,i):
		text = self.ui.saveasSTELLOPT.toPlainText();
		self.fig2.savefig('./'+text, dpi=300)



if __name__ == "__main__":
	app = QApplication(sys.argv) 
	window = MyApp() 
	window.show() 
	sys.exit(app.exec_())
