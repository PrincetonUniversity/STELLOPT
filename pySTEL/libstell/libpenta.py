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
		# NOTE: use endswith (not "in") -- '_call_penta' is also a substring of
		# '_call_penta_surface', so a plain substring match is ambiguous now that
		# both symbols exist in libpenta.so.
		names = [ s for s in attrs if module in s and s.endswith(func)]
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

	def call_PENTA_surface(self, Matom_prof, Zatom_prof,
			ne, dnedrho, te, dtedrho, ni, dnidrho, ti, dtidrho,
			eq_Aminor, eq_Rmajor, vp, chip, phip, iota, btheta, bzeta, bsq,
			DKES_K_surf, rho_k_surf,
			DKES_NUSTAR, DKES_ERSTAR, DKES_D11_surf, DKES_D31_surf, DKES_D33_surf,
			Er_min_Vcm, Er_max_Vcm, EparB_surf, Er_k_surf, Er_root_type, look_for_ambipolar):
		"""Wrapper to call_PENTA_surface in penta_interface_mod

		Computes the ambipolar Er root and neoclassical transport coefficients
		on a SINGLE DKES surface -- i.e. the body of call_PENTA's internal
		"DO k=1,ns_dkes" loop, exposed as its own entry point. This lets
		independent surfaces be computed by independent OS processes (e.g. via
		concurrent.futures.ProcessPoolExecutor), since PENTA keeps its working
		state in Fortran module-level (SAVE) variables that are NOT safe to
		share across concurrent calls within a single process.

		Unlike call_PENTA, this does NOT interpolate onto an output_rho grid --
		callers computing multiple surfaces (e.g. in parallel workers) should
		gather the raw per-surface outputs and call call_PENTA_interpolate once,
		serially, in the parent process.

		Parameters
		----------
		Matom_prof : real (nion_prof)
			Ion species atomic masses [amu].
		Zatom_prof : int (nion_prof)
			Ion species charge numbers [-].
		ne, dnedrho, te, dtedrho : real
			Electron density [m^-3], d(ne)/drho, temperature [eV], d(te)/drho at this surface.
		ni, dnidrho, ti, dtidrho : real (nion_prof)
			Ion density [m^-3], d(ni)/drho, temperature [eV], d(ti)/drho at this surface.
		eq_Aminor, eq_Rmajor : real
			Equilibrium minor/major radius [m].
		vp, chip, phip, iota, btheta, bzeta, bsq : real
			Equilibrium quantities on this DKES surface (VMEC normalization).
		DKES_K_surf : int
			DKES surface index (radial grid point index) of this surface.
		rho_k_surf : real
			Normalized effective radius (r/a) of this DKES surface.
		DKES_NUSTAR : real (ncstar)
			Collisionality grid used by the DKES coefficients.
		DKES_ERSTAR : real (nestar)
			Normalized Er grid used by the DKES coefficients.
		DKES_D11_surf, DKES_D31_surf, DKES_D33_surf : real (ncstar,nestar)
			DKES transport coefficients on this surface.
		Er_min_Vcm, Er_max_Vcm : real
			Search range for the ambipolar Er root [V/cm].
		EparB_surf : real
			<E.B> on this surface.
		Er_k_surf : real
			Er value to use on this surface when look_for_ambipolar is False [V/cm].
		Er_root_type : str
			Which ambipolar root to pick ('ion_root', 'electron_root', 'unstable_root').
		look_for_ambipolar : bool
			If True, search for the ambipolar Er root; if False, use Er_k_surf directly.

		Returns
		-------
		JBS_surf : real
			Bootstrap current density on this surface.
		Er_surf : real
			Ambipolar Er [V/cm] on this surface.
		Dn_surf, cn_surf, Dp_surf, cp_surf : real (nion_prof+1)
			Particle/heat diffusion and convection coefficients on this surface,
			one entry per species (electrons first, then ions in Zatom_prof order).
		"""
		import ctypes as ct
		import numpy as np

		module_name = self.s1+'penta_interface_mod_'+self.s2
		call_penta_surface_h = getattr(self.libpenta, module_name+'_call_penta_surface'+self.s3)

		nion_prof = len(Matom_prof)
		ncstar    = len(DKES_NUSTAR)
		nestar    = len(DKES_ERSTAR)

		call_penta_surface_h.argtypes = [
			ct.POINTER(ct.c_int), ct.POINTER(ct.c_int), ct.POINTER(ct.c_int), \
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
			ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), \
			ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), \
			ct.c_long
		]
		call_penta_surface_h.restype = None

		# Scalars
		nion_prof_c   = ct.c_int(nion_prof)
		ncstar_c      = ct.c_int(ncstar)
		nestar_c      = ct.c_int(nestar)
		ne_c          = ct.c_double(ne)
		dnedrho_c     = ct.c_double(dnedrho)
		te_c          = ct.c_double(te)
		dtedrho_c     = ct.c_double(dtedrho)
		eq_Aminor_c   = ct.c_double(eq_Aminor)
		eq_Rmajor_c   = ct.c_double(eq_Rmajor)
		vp_c          = ct.c_double(vp)
		chip_c        = ct.c_double(chip)
		phip_c        = ct.c_double(phip)
		iota_c        = ct.c_double(iota)
		btheta_c      = ct.c_double(btheta)
		bzeta_c       = ct.c_double(bzeta)
		bsq_c         = ct.c_double(bsq)
		DKES_K_surf_c = ct.c_int(DKES_K_surf)
		rho_k_surf_c  = ct.c_double(rho_k_surf)
		Er_min_Vcm_c  = ct.c_double(Er_min_Vcm)
		Er_max_Vcm_c  = ct.c_double(Er_max_Vcm)
		EparB_surf_c  = ct.c_double(EparB_surf)
		Er_k_surf_c   = ct.c_double(Er_k_surf)
		look_for_ambipolar_c = ct.c_int(1 if look_for_ambipolar else 0)

		# 1D real/int arrays (nion_prof / ncstar / nestar sized)
		Matom_prof  = np.ascontiguousarray(Matom_prof,  dtype=np.float64)
		Zatom_prof  = np.ascontiguousarray(Zatom_prof,  dtype=np.int32)
		ni          = np.ascontiguousarray(ni,          dtype=np.float64)
		dnidrho     = np.ascontiguousarray(dnidrho,     dtype=np.float64)
		ti          = np.ascontiguousarray(ti,          dtype=np.float64)
		dtidrho     = np.ascontiguousarray(dtidrho,     dtype=np.float64)
		DKES_NUSTAR = np.ascontiguousarray(DKES_NUSTAR, dtype=np.float64)
		DKES_ERSTAR = np.ascontiguousarray(DKES_ERSTAR, dtype=np.float64)

		# 2D real arrays for this single surface: Fortran (column-major) order
		DKES_D11_surf = np.asfortranarray(DKES_D11_surf, dtype=np.float64)
		DKES_D31_surf = np.asfortranarray(DKES_D31_surf, dtype=np.float64)
		DKES_D33_surf = np.asfortranarray(DKES_D33_surf, dtype=np.float64)

		# Character(Len=100) dummy argument: pad/truncate to exactly 100 chars.
		Er_root_type_c = Er_root_type.ljust(100)[:100].encode('UTF-8')

		# Outputs
		JBS_surf_c = ct.c_double(0.0)
		Er_surf_c  = ct.c_double(0.0)
		Dn_surf = np.zeros(nion_prof+1, dtype=np.float64)
		cn_surf = np.zeros(nion_prof+1, dtype=np.float64)
		Dp_surf = np.zeros(nion_prof+1, dtype=np.float64)
		cp_surf = np.zeros(nion_prof+1, dtype=np.float64)

		call_penta_surface_h(
			ct.byref(nion_prof_c), ct.byref(ncstar_c), ct.byref(nestar_c), \
			Matom_prof.ctypes.data_as(ct.POINTER(ct.c_double)), \
			Zatom_prof.ctypes.data_as(ct.POINTER(ct.c_int)), \
			ct.byref(ne_c), ct.byref(dnedrho_c), ct.byref(te_c), ct.byref(dtedrho_c), \
			ni.ctypes.data_as(ct.POINTER(ct.c_double)), \
			dnidrho.ctypes.data_as(ct.POINTER(ct.c_double)), \
			ti.ctypes.data_as(ct.POINTER(ct.c_double)), \
			dtidrho.ctypes.data_as(ct.POINTER(ct.c_double)), \
			ct.byref(eq_Aminor_c), ct.byref(eq_Rmajor_c), \
			ct.byref(vp_c), ct.byref(chip_c), ct.byref(phip_c), ct.byref(iota_c), \
			ct.byref(btheta_c), ct.byref(bzeta_c), ct.byref(bsq_c), \
			ct.byref(DKES_K_surf_c), ct.byref(rho_k_surf_c), \
			DKES_NUSTAR.ctypes.data_as(ct.POINTER(ct.c_double)), \
			DKES_ERSTAR.ctypes.data_as(ct.POINTER(ct.c_double)), \
			DKES_D11_surf.ctypes.data_as(ct.POINTER(ct.c_double)), \
			DKES_D31_surf.ctypes.data_as(ct.POINTER(ct.c_double)), \
			DKES_D33_surf.ctypes.data_as(ct.POINTER(ct.c_double)), \
			ct.byref(Er_min_Vcm_c), ct.byref(Er_max_Vcm_c), \
			ct.byref(EparB_surf_c), ct.byref(Er_k_surf_c), \
			Er_root_type_c, ct.byref(look_for_ambipolar_c), \
			ct.byref(JBS_surf_c), ct.byref(Er_surf_c), \
			Dn_surf.ctypes.data_as(ct.POINTER(ct.c_double)), \
			cn_surf.ctypes.data_as(ct.POINTER(ct.c_double)), \
			Dp_surf.ctypes.data_as(ct.POINTER(ct.c_double)), \
			cp_surf.ctypes.data_as(ct.POINTER(ct.c_double)), \
			ct.c_long(len(Er_root_type_c))
		)

		return JBS_surf_c.value, Er_surf_c.value, Dn_surf, cn_surf, Dp_surf, cp_surf

	def call_PENTA_interpolate(self, rho_k, Er_PENTA, Dn_PENTA, cn_PENTA, Dp_PENTA, cp_PENTA, output_rho):
		"""Wrapper to call_PENTA_interpolate in penta_interface_mod

		Interpolates the per-surface ambipolar root and neoclassical transport
		coefficients (as returned by call_PENTA_surface, one call per DKES
		surface, gathered here into arrays indexed by surface) from the DKES
		radial grid (rho_k) onto output_rho. This is the second half of what
		call_PENTA does internally; call it once, serially, after gathering
		every surface's call_PENTA_surface result (e.g. from parallel workers).

		Parameters
		----------
		rho_k : real (ns_dkes)
			Normalized effective radius (r/a) of each DKES surface.
		Er_PENTA : real (ns_dkes)
			Ambipolar Er [V/cm] on each DKES surface.
		Dn_PENTA, cn_PENTA, Dp_PENTA, cp_PENTA : real (nion_prof+1,ns_dkes)
			Particle/heat diffusion and convection coefficients on each DKES
			surface, one row per species (electrons first, then ions).
		output_rho : real (output_nrho)
			Radial grid (r/a) on which to return the transport coefficients.

		Returns
		-------
		output_Er : real (output_nrho)
			Ambipolar Er [V/cm] on output_rho.
		output_Dn, output_cn, output_Dp, output_cp : real (nion_prof+1,output_nrho)
			Particle/heat diffusion and convection coefficients on output_rho.
		"""
		import ctypes as ct
		import numpy as np

		module_name = self.s1+'penta_interface_mod_'+self.s2
		call_penta_interp_h = getattr(self.libpenta, module_name+'_call_penta_interpolate'+self.s3)

		ns_dkes     = len(rho_k)
		nion_prof   = Dn_PENTA.shape[0] - 1
		output_nrho = len(output_rho)

		call_penta_interp_h.argtypes = [
			ct.POINTER(ct.c_int), ct.POINTER(ct.c_int), ct.POINTER(ct.c_int), \
			ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), \
			ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), \
			ct.POINTER(ct.c_double), \
			ct.POINTER(ct.c_double), \
			ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double)
		]
		call_penta_interp_h.restype = None

		ns_dkes_c     = ct.c_int(ns_dkes)
		nion_prof_c   = ct.c_int(nion_prof)
		output_nrho_c = ct.c_int(output_nrho)

		rho_k       = np.ascontiguousarray(rho_k,       dtype=np.float64)
		Er_PENTA    = np.ascontiguousarray(Er_PENTA,    dtype=np.float64)
		output_rho  = np.ascontiguousarray(output_rho,  dtype=np.float64)

		# 2D real arrays: Fortran (column-major) order
		Dn_PENTA = np.asfortranarray(Dn_PENTA, dtype=np.float64)
		cn_PENTA = np.asfortranarray(cn_PENTA, dtype=np.float64)
		Dp_PENTA = np.asfortranarray(Dp_PENTA, dtype=np.float64)
		cp_PENTA = np.asfortranarray(cp_PENTA, dtype=np.float64)

		output_Er = np.zeros(output_nrho, dtype=np.float64)
		output_Dn = np.zeros((nion_prof+1, output_nrho), order='F', dtype=np.float64)
		output_cn = np.zeros((nion_prof+1, output_nrho), order='F', dtype=np.float64)
		output_Dp = np.zeros((nion_prof+1, output_nrho), order='F', dtype=np.float64)
		output_cp = np.zeros((nion_prof+1, output_nrho), order='F', dtype=np.float64)

		call_penta_interp_h(
			ct.byref(ns_dkes_c), ct.byref(nion_prof_c), ct.byref(output_nrho_c), \
			rho_k.ctypes.data_as(ct.POINTER(ct.c_double)), \
			Er_PENTA.ctypes.data_as(ct.POINTER(ct.c_double)), \
			Dn_PENTA.ctypes.data_as(ct.POINTER(ct.c_double)), \
			cn_PENTA.ctypes.data_as(ct.POINTER(ct.c_double)), \
			Dp_PENTA.ctypes.data_as(ct.POINTER(ct.c_double)), \
			cp_PENTA.ctypes.data_as(ct.POINTER(ct.c_double)), \
			output_rho.ctypes.data_as(ct.POINTER(ct.c_double)), \
			output_Er.ctypes.data_as(ct.POINTER(ct.c_double)), \
			output_Dn.ctypes.data_as(ct.POINTER(ct.c_double)), \
			output_cn.ctypes.data_as(ct.POINTER(ct.c_double)), \
			output_Dp.ctypes.data_as(ct.POINTER(ct.c_double)), \
			output_cp.ctypes.data_as(ct.POINTER(ct.c_double))
		)

		return output_Er, output_Dn, output_cn, output_Dp, output_cp


# --- Helpers for running call_PENTA_surface across a persistent
# concurrent.futures.ProcessPoolExecutor. These must be module-level (not
# bound methods) since the ctypes.CDLL handle inside a LIBPENTA instance is
# not picklable -- each worker process instead builds its own LIBPENTA once,
# via _init_NEO_worker as the pool's `initializer`, and keeps it in a
# process-local global.

_worker_libpenta = None

def _init_NEO_worker():
	global _worker_libpenta
	_worker_libpenta = LIBPENTA()

def _call_PENTA_surface_worker(args):
	global _worker_libpenta
	if _worker_libpenta is None:
		# Fallback in case the pool was created without _init_NEO_worker as
		# its initializer.
		_worker_libpenta = LIBPENTA()
	return _worker_libpenta.call_PENTA_surface(*args)
