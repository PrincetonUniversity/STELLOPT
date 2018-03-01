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
                              set_module_var, safe_close
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
		self.indata={}
		self.indata['am']=np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
		# Setup Components
		self.ui.TableArrays.setRowCount(1)
		self.ui.TableArrays.setColumnCount(20)
		self.ui.TableArrays.show()
		# Setup Plot
		self.fig = Figure(figsize=(2,2),dpi=100)
		self.ax = self.fig.add_subplot(111)
		self.canvas = FigureCanvas(self.fig)
		self.ui.Plotbox.addWidget(self.canvas)
		# Callbacks
		self.ui.ButtonLoadIndata.clicked.connect(self.LoadIndata)
		self.ui.ComboBoxArrays.currentIndexChanged.connect(self.UpdateArrays)
		self.ui.ComboBoxPType.currentIndexChanged.connect(self.UpdatePType)
		self.ui.TableArrays.cellChanged.connect(self.DataArrays)

	def LoadIndata(self,i):
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
			#print(num,item)
			self.ui.TableNsArray.setItem(0,num, QTableWidgetItem(str(item)))
			self.ui.TableNsArray.setItem(1,num, QTableWidgetItem(str(niter_array[num])))
			self.ui.TableNsArray.setItem(2,num, QTableWidgetItem(str(ftol_array[num])))
		self.ui.TableNsArray.show()
		self.UpdateArrays()

	def UpdatePType(self):
		data_name=self.ui.ComboBoxArrays.currentText()
		data_name = data_name.lower()
		type_name=self.ui.ComboBoxPType.currentText()
		type_name = type_name.lower()
		if data_name == 'am' or data_name == 'am_aux':
			set_module_var('vmec_input','pmass_type',type_name)
		elif data_name == 'ac' or data_name == 'ac_aux':
			set_module_var('vmec_input','pcurr_type',type_name)
		elif data_name == 'ai' or data_name == 'ai_aux':
			set_module_var('vmec_input','piota_type',type_name)

	def UpdateArrays(self):
		# Updates values in array box
		dex=self.ui.ComboBoxArrays.currentIndex()
		data_name=self.ui.ComboBoxArrays.currentText()
		data_name = data_name.lower()
		# Update PType accordingly
		strtmp = 'none'
		if data_name == 'am' or data_name == 'am_aux':
			strtmp = str(self.indata['pmass_type']).strip()
		if data_name == 'ai' or data_name == 'ai_aux':
			strtmp = str(self.indata['piota_type']).strip()
		if data_name == 'ac' or data_name == 'ac_aux':
			strtmp = str(self.indata['pcurr_type']).strip()
		if strtmp != 'none':
			self.ui.ComboBoxPType.setCurrentIndex(self.ui.ComboBoxPType.findText(strtmp))
		# Update the table
		self.ui.TableArrays.setRowCount(1)
		self.ui.TableArrays.setColumnCount(20)
		if data_name == 'am' or data_name == 'ac' or data_name == 'ai':
			for num,item in enumerate(self.indata[data_name], start=0):
				self.ui.TableArrays.setItem(0,num, QTableWidgetItem(str(item)))
		if data_name == 'am_aux' or data_name == 'ac_aux' or data_name == 'ai_aux':
			aux_mask = self.indata[data_name+'_s']>=0
			self.ui.TableArrays.setRowCount(2)
			if (len(aux_mask) > 0):
				self.ui.TableArrays.setColumnCount(len(aux_mask))
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
		row = self.ui.TableArrays.currentColumn()
		if len(data_name) == 2:
			self.indata[data_name][col]=value
			set_module_var('vmec_input',data_name,self.indata[data_name])
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
		for i in range(99): s[i]=(i-1.0)/(98)
		if data_name[0:2] == 'am':
			for i,xx in enumerate(s):
				f[i] = pmass(xx)
				self.ax.set_ylabel('Pressure')
		elif data_name[0:2] == 'ac':
			for i,xx in enumerate(s):
				f[i] = pcurr(xx)
				self.ax.set_ylabel('Current')
		elif data_name[0:2] == 'ai':
			for i,xx in enumerate(s):
				f[i] = piota(xx)
				self.ax.set_ylabel('Iota')
		self.ax.plot(s,f)
		self.ax.set_xlabel('Norm. Tor. Flux (s)')
		self.ax.set_aspect('auto')
		self.canvas.draw()





if __name__ == "__main__":
	app = QApplication(sys.argv) 
	window = MyApp() 
	window.show() 
	sys.exit(app.exec_())
