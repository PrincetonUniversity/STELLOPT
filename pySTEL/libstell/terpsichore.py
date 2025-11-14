#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and handling 
TERPSICHORE stability information.
"""

# Libraries
from libstell.libstell import FourierRep, LIBSTELL

# Constants

# DIAGNO Class
class TERPSICHORE(FourierRep):
	"""Class for working with TERPSICHORE data

	"""
	def __init__(self):
		self.libterp = None

	def initlibterp(self):
		"""Sets up the libterp libraray

		This routine loads the TERPSICHORE shared library libterpsichore.so
		"""
		import os,sys
		import ctypes as ct
		from subprocess import Popen, PIPE
		self.TERPSICHORE_PATH = os.environ["TERPSICHORE_PATH"]
		self.PATH_TO_TERPSICHORE = os.path.join(self.TERPSICHORE_PATH,'libterpsichore.so')
		try:
			self.libterp = ct.cdll.LoadLibrary(self.PATH_TO_TERPSICHORE)
		except:
			print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
			print("!!  Could not load libraray libterpsichore.so     !!")
			print(f"!!  PATH: {self.PATH_TO_TERPSICHORE}    !!")
			print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
			return
		# Figure out uderscoring
		out = Popen(args="nm "+self.PATH_TO_TERPSICHORE, 
			shell=True, 
			stdout=PIPE).communicate()[0].decode("utf-8")
		attrs = [ i.split(" ")[-1].replace("\r", "") \
			for i in out.split("\n") if " T " in i]
		func = '_read_fort_23'
		module = 'read_terpsichore_mod_' 
		names = [ s for s in attrs if module in s and func in s]
		name = names[0].replace(module, ',')
		name = name.replace(func, ',')
		self.s1, self.s2, self.s3 = name.split(',')
		# Weird OSX behavior
		if self.s1=='___':
			self.s1='__'


	def write_eq_input(self,vmec):
		"""Computes TERPSICHORE input from VMEC data

		This routine takes a VMEC output data class and toroidal mode
		number and computes TERPSICHORE input data.

		Parameters
		----------
		vmec : VMEC Class
			VMEC class containting wout information
		"""
		import numpy as np
		rmu0 = np.pi*4E-7
		# Create the equilibrium file
		f = open('terpsichore_eq.'+vmec.input_extension.strip(),'w')
		f.write(f' {vmec.wb:22.12E} {vmec.gamma:12.5E}{1.0/float(vmec.nfp):12.5E}{0.0:12.5E}')
		f.write(f'{vmec.mnmax_nyq:4d}{vmec.ns:4d}{round(np.max(vmec.xm_nyq)+1):4d}')
		f.write(f'{round(np.max(vmec.xn_nyq)/vmec.nfp):4d}{1:4d}{1:4d}')
		f.write(f'{round(vmec.itfsq):4d}{round(vmec.niter):4d}{0:1d}\n')
		for k in range(vmec.ns):
			for mn in range(vmec.mnmax_nyq):
				sigc=0; tauc=0; pbpc=0; pppc=0; rmnc = 0.0; zmns = 0.0
				if (vmec.xn_nyq[mn] == 0) and (vmec.xm_nyq[mn] == 0):
					sigc = 1; tauc=1
				for mn2 in range(vmec.mnmax):
					if (vmec.xm[mn2] == vmec.xm_nyq[mn]) and (vmec.xn[mn2] == vmec.xn_nyq[mn]):
						rmnc = vmec.rmnc[k,mn2]
						zmns = vmec.zmns[k,mn2]
				f.write(f' {vmec.xm_nyq[mn,0]:22.14E}{-vmec.xn_nyq[mn,0]:22.14E}{rmnc:22.14E}{zmns:22.14E}{vmec.gmnc[k,mn]:22.14E}\n')
				f.write(f' {sigc:22.14E} {tauc:22.14E}{pbpc:22.14E}{pppc:22.14E}\n')
		for k in range(vmec.ns):
			f.write(f' {vmec.iotas[k,0]:22.14E}{vmec.mass[k,0]:22.14E}{rmu0*vmec.pres[k,0]:22.14E}{-vmec.phip[k,0]:22.14E}{vmec.vp[k,0]:22.14E}\n')
		f.close()

	def create_input(self,vmec,n=0,npertmax=3):
		"""Computes TERPSICHORE input from VMEC data

		This routine takes a VMEC output data class and toroidal mode
		number and computes TERPSICHORE input data.

		Parameters
		----------
		vmec : VMEC Class
			VMEC class containting wout information
		n : integer (optional)
			Toroidal mode number (default = 1.0)
		npertmax : integer (optional)
			Number of NFPs to go left and right of n (default = 3)
		"""
		import numpy as np
		rmu0 = np.pi*4E-7
		# Compute some helpers
		qn = 0.0; mm_max = 36; mms = 55
		ni = vmec.ns-1
		ivac = round(vmec.ns/4)
		# Compute max m and n in Boozer transformation
		mm = int(2**np.ceil(np.log2(2*np.max(vmec.xm))))
		mm = min(mm,mm_max)
		nmax = int(2**np.ceil(np.log2(2*np.max(vmec.xn/vmec.nfp))))
		nmin = -nmax
		# Compute max m and n in Stability calculation
		nsmin = -n - npertmax*vmec.nfp
		nsmax =  n + npertmax*vmec.nfp
		# Check n
		nmin = int(min(nsmin,nmin))
		nmax = int(max(nsmax,nmax))
		# Compute number of realspace points to use
		nj = int(2**np.ceil(np.log2(mm)+2))
		nk = int(2**np.ceil(np.log2(nmax)+2))
		nj = max(nj,64)
		nk = max(nk,32)
		# Create Boozer mode matrix
		lfrz = np.zeros((mm_max+1,nmax-nmin+1),dtype=int)
		for n2 in range(nmin,nmax+1):
			if n2 < 0:
				for m in range(1,mm+1): lfrz[m,n2-nmin] = 1
			else:
				for m in range(mm+1): lfrz[m,n2-nmin] = 1
		mlmnb = np.count_nonzero(lfrz)
		# Create Stability mode matrix
		lfrs = np.zeros((mms+1,nsmax-nsmin+1),dtype=int)
		for j in range(-npertmax,npertmax+1):
			n0 = 1
			ntemp = n+j*vmec.nfp
			if ntemp > 0: n0 = 0
			for m in range(n0,mms+1): lfrs[m,ntemp-nsmin] = 1
			n0 = 1
			ntemp =-n+j*vmec.nfp
			if ntemp > 0: n0 = 0
			for m in range(n0,mms+1): lfrs[m,ntemp-nsmin] = 1
		# temporary fix
		#lfrs[51:,:] = 0
		mlmns = np.count_nonzero(lfrs)
		# First create the namelist data
		# Create an input file
		f = open(f'terpsichore_input_{n:02d}','w')
		f.write(f'               {vmec.input_extension.strip()}\n')
		f.write('C\nC        MM  NMIN  NMAX   MMS NSMIN NSMAX NPROCS INSOL\n')
		f.write(f'   {mm:6d}{nmin:6d}{nmax:6d}{mms:6d}{nsmin:6d}{nsmax:6d}{1:6d}{0:6d}\n')
		f.write('C\nC     TABLE OF FOURIER COEFFIENTS FOR BOOZER COORDINATES\n')
		f.write('C     EQUILIBRIUM SETTINGS ARE COMPUTED FROM FIT/VMEC\nC\n')
		f.write('C M=')
		for m in range(37): f.write(f' {np.remainder(m,10):1d}')
		f.write('  N\n')
		for n2 in range(nmin,nmax+1):
			f.write('    ')
			for m in range(mm_max+1):
				f.write(f'{lfrz[m,n2-nmin]:2d}')
			f.write(f' {n2:1d}\n')
		f.write('C\n      LLAMPR      LVMTPR      LMETPR      LFOUPR\n')
		f.write('           0           0           0           0\n')
		f.write('      LLHSPR      LRHSPR      LEIGPR      LEFCPR\n')
		f.write('           9           9           1           1\n')
		f.write('      LXYZPR      LIOTPL      LDW2PL      LEFCPL\n')
		f.write('           0           1           1           1\n')
		f.write('      LCURRF      LMESHP      LMESHV      LITERS\n')
		f.write('           1           1           3           1\n')
		f.write('      LXYZPL      LEFPLS      LEQVPL      LPRESS\n')
		f.write('           1           1           0           2\n')
		f.write('C\nC    PVAC        PARFAC      QONAX        QN         DSVAC       QVAC    NOWALL\n')
		f.write(f'{1.0001:12.4E}{0.00:12.4E}{1./vmec.iotaf[0,0]:12.4E}{qn:12.4E}{1.00:12.4E}{1.0001:12.4E}     {-2:2d}\n')
		f.write('C\nC    AWALL       EWALL       DWALL       GWALL       DRWAL       DZWAL   NPWALL\n')
		f.write(f'{2.00:12.4E}{1.00:12.4E}{0.50:12.4E}{vmec.rmnc[0,vmec.mn00]:12.4E}{0.00:12.4E}{0.00:12.4E}     {vmec.nfp:2d}\n')
		f.write('C\nC    RPLMIN       XPLO      DELTAJP       WCT      CURFAC\n')
		f.write(f'{1E-5:12.4E}{1E-6:12.4E}{4E-2:12.4E}{vmec.rmnc[0,1]:12.4E}{1.00:12.4E}\n')
		f.write(f'C\nC                                                             MODELK = {1:6d}\n')
		f.write(f'C\nC     NUMBER OF EQUILIBRIUM FIELD PERIODS PER STABILITY PERIOD: NSTA = {vmec.nfp:6d}\n')
		f.write('C\nC     TABLE OF FOURIER COEFFIENTS FOR STABILITY DISPLACEMENTS\n')
		f.write('C\nC M=')
		for m in range(mms+1): f.write(f' {np.remainder(m,10):1d}')
		f.write('  N\n')
		for n2 in range(nsmin,nsmax+1):
			f.write('    ')
			for m in range(mms+1):
				f.write(f'{lfrs[m,n2-nsmin]:2d}')
			f.write(f' {n2:1d}\n')
		f.write('C\nC   NEV NITMAX         AL0     EPSPAM IGREEN MPINIT\n')
		f.write(f'{1:7d}{1500:7d}{-5E-3:12.3E}{1E-4:12.3E}{0:7d}{0:7d}\nC\n')
		f.close()
		# Comput LSSL
		mmaxdf=2*mm
		nmaxdf=2*max(abs(nmin),nmax)
		lssl = self._compute_lssl(mm,nmin,nmax,mmaxdf,nmaxdf,lfrz)
		# Compute LSSD 
		mmaxdf=2*mms
		nmaxdf=2*max(abs(nsmin),nsmax)
		lssd = self._compute_lssd(mms,nsmin,nsmax,mmaxdf,nmaxdf,lfrs,vmec.nfp)
		# Now get mmax and nmax
		mmaxdf=max(2*mms,2*mm)
		nmaxdf=max(2*max(abs(nsmin),nsmax),2*max(abs(nmin),nmax))
		# Output information for module file
		print(f'!      {vmec.input_extension.strip()} (n={n:2d})')
		print(f'       INTEGER :: NI = {ni:4d}')
		print(f'       INTEGER :: IVAC = {ivac:4d}')
		print(f'       INTEGER :: NVI = {ni+ivac:4d}')
		print(f'       INTEGER :: NJ = {nj:4d}')
		print(f'       INTEGER :: NK = {nk:4d}')
		print(f'       INTEGER :: NJK = {nj*nk:6d}')
		print(f'       INTEGER :: MLMNV = {vmec.mnmax_nyq:4d}')
		print(f'       INTEGER :: MLMNB = {mlmnb:4d}')
		print(f'       INTEGER :: LSSL = {lssl:4d}')
		print(f'       INTEGER :: MMAXDF = {mmaxdf:6d}')
		print(f'       INTEGER :: NMAXDF = {nmaxdf:6d}')
		print(f'       INTEGER :: ND = {ni+ivac:4d}')
		print(f'       INTEGER :: ND1 = {ni+ivac+1:4d}')
		print(f'       INTEGER :: LSSD = {lssd:4d}')
		print(f'       INTEGER :: MLMNS = {mlmns:4d}')
		print(f'       INTEGER :: MD = {mlmns:4d}')
		print(f'       INTEGER :: MDY = {mlmns:4d}')
		print(f'       INTEGER :: NA = {2*mlmns*(ni+ivac)+mlmns:6d}')

	def _compute_lssl(self,mm,nmin,nmax,mmaxdf,nmaxdf,lfrz):
		"""Computes the LSSL

		Computes the LSSL term tprgl0.module_ap.f line 883
		"""
		import numpy as np
		# Compute ML and NL
		lmnl = 0; ml=[]; nl=[]
		for m in range(mm+1):
			for n in range(nmin,nmax+1):
				if lfrz[m,n-nmin]>0:
					if (m != 0) or (n != 0):
						lmnl = lmnl + 1
						ml.append(m)
						nl.append(n)
		# Now compute lss mxdif=[-M,M] nxdif=[-2N,2N]
		lfx = np.zeros((2*mmaxdf+1,4*nmaxdf+1),dtype=int)
		lss = 0
		for lc in range(lmnl):
			for lr in range(lmnl):
				mxdif = ml[lc] - ml[lr]
				nxdif = nl[lc] - nl[lr]
				mdex = mxdif + mmaxdf
				ndex = nxdif + 2*nmaxdf
				if (lfx[mdex,ndex] <= 0):
					lfx[mdex,ndex] = 1
					lss = lss + 1
		for lc in range(lmnl):
			for lr in range(lmnl):
				mxdif = ml[lc] + ml[lr]
				nxdif = nl[lc] + nl[lr]
				mdex = mxdif + mmaxdf
				ndex = nxdif + 2*nmaxdf
				if (lfx[mdex,ndex] <= 0):
					lfx[mdex,ndex] = 1
					lss = lss + 1
		return lss

	def _compute_lssd(self,mm,nmin,nmax,mmaxdf,nmaxdf,lfrz,nfp):
		"""Computes the LSSD

		Computes the LSSD term tprgl0.module_ap.f line 1412
		"""
		import numpy as np
		# Compute ML and NL
		lmns = 0; ms=[]; ns=[]
		for n in range(nmin,nmax+1):
			for m in range(mm+1):
				if lfrz[m,n-nmin]>0:
					lmns = lmns + 1
					ms.append(m)
					ns.append(n)
		# Now compute lss mxdif=[-M,M] nxdif=[-2N,2N]
		lfx = np.zeros((2*mmaxdf+1,4*nmaxdf+1),dtype=int)
		lss = 0
		for lc in range(lmns):
			for lr in range(lmns):
				mxdif = ms[lc] - ms[lr]
				nxdif = ns[lc] - ns[lr]
				mdex = mxdif + mmaxdf
				ndex = nxdif + 2*nmaxdf
				if (np.mod(nxdif,nfp) == 0):
					if (lfx[mdex,ndex] <= 0):
						lfx[mdex,ndex] = 1
						lss = lss + 1
		for lc in range(lmns):
			for lr in range(lmns):
				mxdif = ms[lc] + ms[lr]
				nxdif = ns[lc] + ns[lr]
				mdex = mxdif + mmaxdf
				ndex = nxdif + 2*nmaxdf
				if (np.mod(nxdif,nfp) == 0):
					if (lfx[mdex,ndex] <= 0):
						lfx[mdex,ndex] = 1
						lss = lss + 1
		return lss

	def read_terpsichore_17(self,filename='fort.17'):
		"""Reads the TERPSICHORE fort.17 file.

		This routine reads the TERPSICHORE fort.17 file containing
		the equilibrium data.

		Parameters
		----------
		filename : str (optional)
			File to read (default: fort.17)
		"""
		import numpy as np
		# Load file
		f = open(filename,'r')
		lines = f.readlines()
		f.close()
		# Get Input Header
		head_txt = lines[3].split()
		lmap = int(head_txt[2])
		ni   = int(head_txt[5])+1
		lmn  = int(head_txt[7][1:])
		nimn = ni*lmn
		nimn6 = nimn+6
		# Read Input harmonics
		xm = []; xn=[]; rmnc=[]; zmns=[]
		for line in lines[6:nimn6]:
			(i_txt,l_txt,m_txt,n_txt,r_txt,z_txt) = line.split()
			xm.append(int(m_txt))
			xn.append(int(n_txt))
			rmnc.append(float(r_txt))
			zmns.append(float(z_txt))
		self.xm_in = np.array([xm[0:lmn]]).T
		self.xn_in = -np.array([xn[0:lmn]]).T # Need in (mu+nv) format
		self.rmnc_in = np.reshape(rmnc,(ni,lmn))
		self.zmns_in = np.reshape(zmns,(ni,lmn))
		temp = abs(self.xn_in)
		self.nfp_in = np.min(temp[temp>0])
		# Read transformed quantities
		if len(lines) > nimn6:
			nimn66 = nimn6+6
			head_txt = lines[nimn6+3].split()
			#TBD copy from above
			lmap = int(head_txt[2])
			#ni   = int(head_txt[5])+1 # should be the same
			lmnb  = int(head_txt[7][1:])
			nimnb = ni*lmnb
			nimnb6 = nimn + nimn66
			# Read TRANSFORMED harmonics
			xm = []; xn=[]; rmnc=[]; zmns=[]
			for line in lines[nimn66:]:
				(i_txt,l_txt,m_txt,n_txt,r_txt,z_txt) = line.split()
				xm.append(int(m_txt))
				xn.append(int(n_txt))
				rmnc.append(float(r_txt))
				zmns.append(float(z_txt))
			self.xm_b = np.array([xm[0:lmnb]]).T
			self.xn_b = np.array([xn[0:lmnb]]).T
			self.rmnc_b = np.reshape(rmnc,(ni,lmnb))
			self.zmns_b = np.reshape(zmns,(ni,lmnb))
			temp = abs(self.xn_b)
			self.nfp_b = np.min(temp[temp>0])

	def read_terpsichore_19(self,filename='fort.19'):
		"""Reads the TERPSICHORE fort.19 file.

		This routine reads the TERPSICHORE fort.19 file containing
		the plasma-vacuum interface information and wall data.

		Parameters
		----------
		filename : str (optional)
			File to read (default: fort.19)
		"""
		import numpy as np
		# Load file
		f = open(filename,'r')
		lines = f.readlines()
		f.close()
		# Get header
		[txt1,txt2,txt3] = lines[0].split()
		self.njk = int(txt1)
		self.nj = int(txt2)
		self.nk = int(txt3)
		# Get data
		rwall = []; zwall = []; rpvi=[]; zpvi=[]
		for line in lines[1:]:
			[txt1,txt2,txt3,txt4,txt5] = line.split()
			rwall.append(float(txt2))
			zwall.append(float(txt3))
			rpvi.append(float(txt4))
			zpvi.append(float(txt5))
		self.Rwall = np.reshape(rwall,(self.nj,self.nk),order='F')
		self.Zwall = np.reshape(zwall,(self.nj,self.nk),order='F')
		self.Rpvi = np.reshape(rpvi,(self.nj,self.nk),order='F')
		self.Zpvi = np.reshape(zpvi,(self.nj,self.nk),order='F')
		# Output doesn't close curves
		self.Rwall = np.vstack((self.Rwall,self.Rwall[0,:]))
		self.Zwall = np.vstack((self.Zwall,self.Zwall[0,:]))
		self.Rpvi = np.vstack((self.Rpvi,self.Rpvi[0,:]))
		self.Zpvi = np.vstack((self.Zpvi,self.Zpvi[0,:]))

	def read_terpsichore_22(self,filename='fort.22'):
		"""Reads the TERPSICHORE fort.22 file.

		This routine reads the TERPSICHORE fort.22 

		Parameters
		----------
		filename : str (optional)
			File to read (default: fort.22)
		"""
		import numpy as np
		# Load file
		f = open(filename,'r')
		lines = f.readlines()
		f.close()
		# Get header
		[txt1,txt2,txt3] = lines[0].split()
		self.ni = int(txt1)
		self.lmns = int(txt2)
		self.gamma = float(txt3)
		# Get pth
		f = []; nline = 1
		while (len(f) < self.ni):
			txt = lines[nline].split()
			for t in txt:
				f.append(float(t))
			nline = nline + 1
		self.pth = np.array(f)
		# Get aiota
		f = []
		while (len(f) < self.ni):
			txt = lines[nline].split()
			for t in txt:
				f.append(float(t))
			nline = nline + 1
		self.aiota = np.array(f)
		# Get wpsi
		f = []
		while (len(f) < self.ni):
			txt = lines[nline].split()
			for t in txt:
				f.append(float(t))
			nline = nline + 1
		self.wpsi = np.array(f)
		# Get n
		f = []
		while (len(f) < self.lmns*2):
			txt = lines[nline].split()
			for t in txt:
				f.append(int(t))
			nline = nline + 1
		self.ms=f[0:self.lmns]
		self.ns=f[self.lmns:]
		# Loop over am,pvp,pvpi
		am=[]; pvp=[]; pvpi=[]
		for ltemp in range(self.lmns):
			# Get am
			f = []
			while (len(f) < self.ni):
				txt = lines[nline].split()
				for t in txt:
					f.append(float(t))
				nline = nline + 1
			am.append(f)
			# Get pvp
			f = []
			while (len(f) < self.ni):
				txt = lines[nline].split()
				for t in txt:
					f.append(float(t))
				nline = nline + 1
			pvp.append(f)
			# Get pvpi
			f = []
			while (len(f) < self.ni):
				txt = lines[nline].split()
				for t in txt:
					f.append(float(t))
				nline = nline + 1
			pvpi.append(f)
		self.am  = np.reshape(am,(self.lmns,self.ni))
		self.pvp  = np.reshape(pvp,(self.lmns,self.ni))
		self.pvpi  = np.reshape(pvpi,(self.lmns,self.ni))

	def read_terpsichore_23(self,filename='fort.23'):
		"""Reads the TERPSICHORE fort.23 file.

		This routine reads the TERPSICHORE fort.23

		Parameters
		----------
		filename : str (optional)
			File to read (default: fort.23)
		"""
		import ctypes as ct
		# Attempt to load the libraray
		if type(self.libterp) == type(None):
			self.initlibterp()
		if type(self.libterp) == type(None): return
		# We use an added routine as a helper
		module_name = self.s1+'read_terpsichore_mod_'+self.s2
		read_fort_23 = getattr(self.libterp,module_name+'_read_fort_23'+self.s3)
		read_fort_23.argtypes = [ct.c_char_p,ct.POINTER(ct.c_int),ct.c_long]
		read_fort_23.restype=None
		istat = ct.c_int(0)
		read_fort_23(filename.encode('UTF-8'),ct.byref(istat),len(filename))
		if not (istat.value == 0):
			return
		# Get Scalars
		intList  = ['ni','nj','nk','njk','nsta','nper','lmns','nvi','lmnb','modelk']
		intLen=[1]*len(intList)
		realList=['curfac', 'wp', 'wk','parity']
		realLen=[1]*len(realList)
		scalar_data = self._get_module_vars(module_name,intVar=intList,intLen=intLen,realVar=realList,realLen=realLen)
		ni = scalar_data['ni']
		nj = scalar_data['nj']
		nk = scalar_data['nk']
		njk = scalar_data['njk']
		nvi = scalar_data['nvi']
		mlmns = scalar_data['lmns']
		mlmnb = scalar_data['lmnb']
		# Get multi-dimensional data
		intList  = ['ms','ns','mb','nb']
		intLen=[(mlmns,1),(mlmns,1),(mlmnb,1),(mlmnb,1)]
		realList=['s','pth','fpp','ftp','cj','ci','pp', \
			'ql','xi','eta','rmu',\
			'r','z','phv','rs','zs','rt','zt','rp','zp',\
			'sigbs','bjac','gssl','gstl','gttl',\
			'fbjac','bjacs','ftpp','fppp','cip','cjp',\
			'gparp','gperp','sigmab','taub','parkur','parjp']
		realLen=[(nvi+1,1),(ni,1),(ni+1,1),(ni+1,1),(ni,1),(ni,1),(ni,1),\
			(mlmns,1),(nvi+1,mlmns),(nvi,mlmns),(nvi,mlmns),\
			(ni+1,njk),(ni+1,njk),(ni,njk),(ni,njk),(ni,njk),(ni,njk),(ni,njk),(ni,njk),(ni,njk),\
			(ni+1,njk),(nvi+1,njk),(nvi+1,njk),(nvi+1,njk),(nvi+1,njk),\
			(ni,mlmnb),(ni+1,njk),(ni+1,1),(ni+1,1),(ni,1),(ni,1),\
			(ni,njk),(ni,njk),(ni,njk),(ni,njk),(ni,njk),(ni,njk)]
		array_data = self._get_module_vars(module_name,intVar=intList,intLen=intLen,realVar=realList,realLen=realLen)
		# Set the class attributes
		for key in scalar_data:
			setattr(self, key, scalar_data[key])
		for key in array_data:
			setattr(self, key, array_data[key])
		return

	def _get_module_vars(self,modName,booVar=None,booLen=None,\
		intVar=None,intLen=None,realVar=None,realLen=None,\
		charVar=None,charLen=None,\
		ldefined_size_arrays=False):
		"""Reads module variables into a python dictionary

		This rountine helps to streamline reading of module variables
		into a python dictionary. The variable lists (booVar, intVar,
		and realVar) are lists of strings referencing the variable names.
		The length lists (booLen,intLen,realLen) are a list of tupules
		defining the size of the variable.

		Parameters
		----------
		modName : str
			Name of the module eg 'read_wout_mod'
		booVar : list (optional)
			List of strings referencing boolean variables to pull
		booLen : list (optional)
			List of tuples defining array size
		intVar : list (optional)
			List of strings referencing integer variables to pull
		intLen : list (optional)
			List of tuples defining array size
		realVar : list (optional)
			List of strings referencing integer variables to pull
		realLen : list (optional)
			List of tuples defining array size
		ldefined_size_arrays : logical (optional)
			Set to true when reading a predefined size arrays temp(0:20)
		Returns
		-------
		out_data : dict
			Dictionary of module variables
		"""
		import ctypes as ct
		import numpy.ctypeslib as npct
		from math import prod
		out_data={}
		# Booleans
		if booVar:
			ftemp = ct.POINTER(ct.c_bool)
			for i,temp in enumerate(booVar):
				#print(temp,booLen[i])
				if booLen[i]==1:
					out_data[temp]=ct.c_bool.in_dll(self.libterp,modName+'_'+temp+self.s3).value
				else:
					# This works because fortran has 4 byte sized booleans
					if ldefined_size_arrays : ftemp=ct.c_int*prod(booLen[i])
					out_data[temp]=npct.as_array(ftemp.in_dll(self.libterp,modName+'_'+temp+self.s3),booLen[i])>0
		# Integers
		if intVar:
			ftemp = ct.POINTER(ct.c_int)
			for i,temp in enumerate(intVar):
				#print(temp,intLen[i])
				if intLen[i]==1:
					out_data[temp]=ct.c_int.in_dll(self.libterp,modName+'_'+temp+self.s3).value
				else:
					if ldefined_size_arrays : ftemp=ct.c_int*prod(intLen[i])
					out_data[temp]=npct.as_array(ftemp.in_dll(self.libterp,modName+'_'+temp+self.s3),intLen[i])
		# Reals
		if realVar:
			ftemp = ct.POINTER(ct.c_float)
			for i,temp in enumerate(realVar):
				#print(temp,realLen[i])
				if realLen[i]==1:
					out_data[temp]=ct.c_float.in_dll(self.libterp,modName+'_'+temp+self.s3).value
				else:
					if ldefined_size_arrays : ftemp=ct.c_float*prod(realLen[i])
					out_data[temp]=npct.as_array(ftemp.in_dll(self.libterp,modName+'_'+temp+self.s3),realLen[i])
		# Characters
		if charVar:
			ftemp = ct.POINTER(ct.c_char)
			for i,temp in enumerate(charVar):
				#print(temp,charLen[i])
				if charLen[i]==1:
					out_data[temp]=ct.c_char.in_dll(self.libterp,modName+'_'+temp+self.s3).value.decode('UTF-8')
				else:
					if ldefined_size_arrays : ftemp=ct.c_char*prod(charLen[i])
					out_data[temp]=ftemp.in_dll(self.libterp,modName+'_'+temp+self.s3).value.decode('UTF-8')
		return out_data






# Main routine
if __name__=="__main__":
	import sys
	sys.exit(0)