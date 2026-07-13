"""
This library provides a python class for interfacing to libpenta
"""

# Libraries

# Constants

# LIBPENTA Class
class LIBPENTA():
	"""Class for working with PENTA library routines (fortran interfaces via Ctypes)

	"""
	def __init__(self, parent=None):
		import os
		import ctypes as ct
		from subprocess import Popen, PIPE
		self.STELLOPT_PATH = os.environ["STELLOPT_PATH"]
		self.PATH_TO_LIBPENTA = os.path.join(self.STELLOPT_PATH,'PENTA','Release','libpenta.so')
		try:
			self.libpenta = ct.cdll.LoadLibrary(self.PATH_TO_LIBPENTA)
		except:
			print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
			print("!!  Could not load shared libraray libpenta.so    !!")
			print(f"!!  PATH: {self.PATH_TO_LIBPENTA}    !!")
			print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
		# Figure out underscoring
		out = Popen(args="nm "+self.PATH_TO_LIBPENTA,
			shell=True,
			stdout=PIPE).communicate()[0].decode("utf-8")
		attrs = [ i.split(" ")[-1].replace("\r", "") \
			for i in out.split("\n") if " T " in i]
		func = '_call_penta'
		module = 'penta_interface_mod_'
		names = [ s for s in attrs if module in s and func in s]
		name = names[0].replace(module, ',')
		name = name.replace(func, ',')
		self.s1, self.s2, self.s3 = name.split(',')
		# Weird OSX behavior
		if self.s1=='___':
			self.s1='__'

	def call_PENTA(self, Matom_prof, Zatom_prof,
			ne, dnedrho, te, dtedrho, ni, dnidrho, ti, dtidrho,
			eq_Aminor, eq_Rmajor, vp, chip, phip, iota, btheta, bzeta, bsq,
			DKES_K, rho_k,
			DKES_NUSTAR, DKES_ERSTAR, DKES_D11, DKES_D31, DKES_D33,
			Er_min_Vcm, Er_max_Vcm, EparB, Er_k, Er_root_type, look_for_ambipolar,
			output_rho):
		"""Wrapper to call_PENTA in penta_interface_mod

		This routine wrappers the call_PENTA subroutine found in
		PENTA/Sources/penta_interface_mod.f90. It loops PENTA over the
		DKES surfaces described by the per-surface (ns_dkes) inputs,
		finds the ambipolar Er root at each surface and returns the
		bootstrap-related neoclassical transport coefficients.

		Parameters
		----------
		Matom_prof : real (nion_prof)
			Ion species atomic masses [amu].
		Zatom_prof : int (nion_prof)
			Ion species charge numbers [-].
		ne, dnedrho, te, dtedrho : real (ns_dkes)
			Electron density [m^-3], d(ne)/drho, temperature [eV], d(te)/drho.
		ni, dnidrho, ti, dtidrho : real (ns_dkes,nion_prof)
			Ion density [m^-3], d(ni)/drho, temperature [eV], d(ti)/drho.
		eq_Aminor, eq_Rmajor : real
			Equilibrium minor/major radius [m].
		vp, chip, phip, iota, btheta, bzeta, bsq : real (ns_dkes)
			Equilibrium quantities on the DKES surfaces (VMEC normalization).
		DKES_K : int (ns_dkes)
			DKES surface index (radial grid point index) of each entry.
		rho_k : real (ns_dkes)
			Normalized effective radius (r/a) of each DKES surface.
		DKES_NUSTAR : real (ncstar)
			Collisionality grid used by the DKES coefficients.
		DKES_ERSTAR : real (nestar)
			Normalized Er grid used by the DKES coefficients.
		DKES_D11, DKES_D31, DKES_D33 : real (ns_dkes,ncstar,nestar)
			DKES transport coefficients.
		Er_min_Vcm, Er_max_Vcm : real
			Search range for the ambipolar Er root [V/cm].
		EparB : real (ns_dkes)
			<E.B> on each surface.
		Er_k : real (ns_dkes)
			Er value to use on each surface when look_for_ambipolar is False [V/cm].
		Er_root_type : str
			Which ambipolar root to pick ('ion_root', 'electron_root', 'unstable_root').
		look_for_ambipolar : bool
			If True, search for the ambipolar Er root; if False, use Er_k directly.
		output_rho : real (output_nrho)
			Radial grid (r/a) on which to return the transport coefficients.

		Returns
		-------
		output_Er : real (output_nrho)
			Ambipolar Er [V/cm] on output_rho.
		output_Dn, output_cn, output_Dp, output_cp : real (nion_prof+1,output_nrho)
			Particle/heat diffusion and convection coefficients on output_rho,
			one row per species (electrons first, then ions in Zatom_prof order).
		"""
		import ctypes as ct
		import numpy as np

		module_name = self.s1+'penta_interface_mod_'+self.s2
		call_penta_h = getattr(self.libpenta, module_name+'_call_penta'+self.s3)

		# Sizes are inferred from the arrays themselves.
		ns_dkes     = len(ne)
		nion_prof   = len(Matom_prof)
		ncstar      = len(DKES_NUSTAR)
		nestar      = len(DKES_ERSTAR)
		output_nrho = len(output_rho)

		call_penta_h.argtypes = [
			ct.POINTER(ct.c_int), ct.POINTER(ct.c_int), ct.POINTER(ct.c_int), \
			ct.POINTER(ct.c_int), ct.POINTER(ct.c_int), \
			ct.POINTER(ct.c_double), ct.POINTER(ct.c_int), \
			ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), \
			ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), \
			ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), \
			ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), \
			ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), \
			ct.POINTER(ct.c_int), ct.POINTER(ct.c_double), \
			ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), \
			ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), \
			ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), \
			ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), \
			ct.c_char_p, ct.POINTER(ct.c_int), \
			ct.POINTER(ct.c_double), \
			ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), \
			ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), \
			ct.c_long
		]
		call_penta_h.restype = None

		# Scalars
		ns_dkes_c     = ct.c_int(ns_dkes)
		ncstar_c      = ct.c_int(ncstar)
		nestar_c      = ct.c_int(nestar)
		nion_prof_c   = ct.c_int(nion_prof)
		output_nrho_c = ct.c_int(output_nrho)
		eq_Aminor_c   = ct.c_double(eq_Aminor)
		eq_Rmajor_c   = ct.c_double(eq_Rmajor)
		Er_min_Vcm_c  = ct.c_double(Er_min_Vcm)
		Er_max_Vcm_c  = ct.c_double(Er_max_Vcm)
		look_for_ambipolar_c = ct.c_int(1 if look_for_ambipolar else 0)

		# 1D real arrays (ns_dkes / nion_prof / ncstar / nestar / output_nrho sized)
		Matom_prof = np.ascontiguousarray(Matom_prof, dtype=np.float64)
		ne         = np.ascontiguousarray(ne,         dtype=np.float64)
		dnedrho    = np.ascontiguousarray(dnedrho,    dtype=np.float64)
		te         = np.ascontiguousarray(te,         dtype=np.float64)
		dtedrho    = np.ascontiguousarray(dtedrho,    dtype=np.float64)
		vp         = np.ascontiguousarray(vp,         dtype=np.float64)
		chip       = np.ascontiguousarray(chip,       dtype=np.float64)
		phip       = np.ascontiguousarray(phip,       dtype=np.float64)
		iota       = np.ascontiguousarray(iota,       dtype=np.float64)
		btheta     = np.ascontiguousarray(btheta,     dtype=np.float64)
		bzeta      = np.ascontiguousarray(bzeta,      dtype=np.float64)
		bsq        = np.ascontiguousarray(bsq,        dtype=np.float64)
		rho_k      = np.ascontiguousarray(rho_k,      dtype=np.float64)
		DKES_NUSTAR = np.ascontiguousarray(DKES_NUSTAR, dtype=np.float64)
		DKES_ERSTAR = np.ascontiguousarray(DKES_ERSTAR, dtype=np.float64)
		EparB      = np.ascontiguousarray(EparB,      dtype=np.float64)
		Er_k       = np.ascontiguousarray(Er_k,       dtype=np.float64)
		output_rho = np.ascontiguousarray(output_rho, dtype=np.float64)

		# 1D int arrays
		Zatom_prof = np.ascontiguousarray(Zatom_prof, dtype=np.int32)
		DKES_K     = np.ascontiguousarray(DKES_K,     dtype=np.int32)

		# 2D/3D real arrays: Fortran (column-major) order
		ni      = np.asfortranarray(ni,      dtype=np.float64)
		dnidrho = np.asfortranarray(dnidrho, dtype=np.float64)
		ti      = np.asfortranarray(ti,      dtype=np.float64)
		dtidrho = np.asfortranarray(dtidrho, dtype=np.float64)
		DKES_D11 = np.asfortranarray(DKES_D11, dtype=np.float64)
		DKES_D31 = np.asfortranarray(DKES_D31, dtype=np.float64)
		DKES_D33 = np.asfortranarray(DKES_D33, dtype=np.float64)

		# Character(Len=100) dummy argument: pad/truncate to exactly 100 chars.
		Er_root_type_c = Er_root_type.ljust(100)[:100].encode('UTF-8')

		# Outputs (Fortran-ordered so PENTA fills them in the expected layout)
		output_Er = np.zeros(output_nrho, dtype=np.float64)
		output_Dn = np.zeros((nion_prof+1, output_nrho), order='F', dtype=np.float64)
		output_cn = np.zeros((nion_prof+1, output_nrho), order='F', dtype=np.float64)
		output_Dp = np.zeros((nion_prof+1, output_nrho), order='F', dtype=np.float64)
		output_cp = np.zeros((nion_prof+1, output_nrho), order='F', dtype=np.float64)

		call_penta_h(
			ct.byref(ns_dkes_c), ct.byref(ncstar_c), ct.byref(nestar_c), \
			ct.byref(nion_prof_c), ct.byref(output_nrho_c), \
			Matom_prof.ctypes.data_as(ct.POINTER(ct.c_double)), \
			Zatom_prof.ctypes.data_as(ct.POINTER(ct.c_int)), \
			ne.ctypes.data_as(ct.POINTER(ct.c_double)), \
			dnedrho.ctypes.data_as(ct.POINTER(ct.c_double)), \
			te.ctypes.data_as(ct.POINTER(ct.c_double)), \
			dtedrho.ctypes.data_as(ct.POINTER(ct.c_double)), \
			ni.ctypes.data_as(ct.POINTER(ct.c_double)), \
			dnidrho.ctypes.data_as(ct.POINTER(ct.c_double)), \
			ti.ctypes.data_as(ct.POINTER(ct.c_double)), \
			dtidrho.ctypes.data_as(ct.POINTER(ct.c_double)), \
			ct.byref(eq_Aminor_c), ct.byref(eq_Rmajor_c), \
			vp.ctypes.data_as(ct.POINTER(ct.c_double)), \
			chip.ctypes.data_as(ct.POINTER(ct.c_double)), \
			phip.ctypes.data_as(ct.POINTER(ct.c_double)), \
			iota.ctypes.data_as(ct.POINTER(ct.c_double)), \
			btheta.ctypes.data_as(ct.POINTER(ct.c_double)), \
			bzeta.ctypes.data_as(ct.POINTER(ct.c_double)), \
			bsq.ctypes.data_as(ct.POINTER(ct.c_double)), \
			DKES_K.ctypes.data_as(ct.POINTER(ct.c_int)), \
			rho_k.ctypes.data_as(ct.POINTER(ct.c_double)), \
			DKES_NUSTAR.ctypes.data_as(ct.POINTER(ct.c_double)), \
			DKES_ERSTAR.ctypes.data_as(ct.POINTER(ct.c_double)), \
			DKES_D11.ctypes.data_as(ct.POINTER(ct.c_double)), \
			DKES_D31.ctypes.data_as(ct.POINTER(ct.c_double)), \
			DKES_D33.ctypes.data_as(ct.POINTER(ct.c_double)), \
			ct.byref(Er_min_Vcm_c), ct.byref(Er_max_Vcm_c), \
			EparB.ctypes.data_as(ct.POINTER(ct.c_double)), \
			Er_k.ctypes.data_as(ct.POINTER(ct.c_double)), \
			Er_root_type_c, ct.byref(look_for_ambipolar_c), \
			output_rho.ctypes.data_as(ct.POINTER(ct.c_double)), \
			output_Er.ctypes.data_as(ct.POINTER(ct.c_double)), \
			output_Dn.ctypes.data_as(ct.POINTER(ct.c_double)), \
			output_cn.ctypes.data_as(ct.POINTER(ct.c_double)), \
			output_Dp.ctypes.data_as(ct.POINTER(ct.c_double)), \
			output_cp.ctypes.data_as(ct.POINTER(ct.c_double)), \
			ct.c_long(len(Er_root_type_c))
		)

		return output_Er, output_Dn, output_cn, output_Dp, output_cp
