#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for reading and plotting the
convergence-diagnostic files written by the cyl2flx / newt2d inverse
map (R,Z)->(s,u) in LIBSTELL/Sources/Modules/vmec_utils.f when that
module is compiled with -DDEBUG_CYL2FLX.

Each MPI rank writes a file 'cyl2flx_diag_<rank>.dat' with one row per
grid point that failed to converge cleanly (info=-1) and, optionally,
per clear-outside verdict (info=-3, if cyl2flx_log_outside=.TRUE.).

Columns (whitespace separated, header line begins with '#'):
    0  R_target   cylindrical R of the target point          [m]
    1  Nphi       toroidal angle * field periods (nfp*phi)    [rad]
    2  Z_target   cylindrical Z of the target point           [m]
    3  s          solved normalized toroidal flux             [-]
    4  u          solved poloidal angle                       [rad]
    5  fmin       normalized residual (R-Rt)^2+(Z-Zt)^2 * fnorm
    6  dist       geometric matching error sqrt(fmin/fnorm)   [m]
    7  nfe        number of function evaluations used
    8  info       exit code: -1 not converged, -3 outside (s>2)
    9  ndegen     number of degenerate-Jacobian hits (tau~0)
   10  taumin     smallest |Jacobian| Ru*Zs-Zu*Rs encountered [m^2]

The point of the diagnostics is to decide, per equilibrium, whether
failures are driven by an ill-conditioned equilibrium (the map is
near-singular: ndegen>0 / tiny taumin, typically at bean tips and the
X-point) or by a tolerance/iteration-budget problem (nfe pinned at
niter with dist still decreasing, often showing as constant-Z bands).
"""

# Column indices
IR, INPHI, IZ, IS, IU, IFMIN, IDIST, INFE, IINFO, INDEGEN, ITAUMIN = range(11)
COLNAMES = ['R', 'Nphi', 'Z', 's', 'u', 'fmin', 'dist',
            'nfe', 'info', 'ndegen', 'taumin']

# Failure-mode labels and plot colors
MODE_LABELS = ['rescued', 'degenerate', 'starved', 'failed', 'outside']
MODE_COLORS = {'degenerate': 'tab:red',
               'starved':    'tab:orange',
               'rescued':    'tab:green',
               'failed':     'tab:purple',
               'outside':    'tab:gray'}


class Cyl2FlxDiag():
	"""Class for reading and plotting cyl2flx convergence diagnostics

	"""
	def __init__(self):
		self.data = None    # ndarray (N, 11), one row per logged point
		self.files = []     # list of files that were read

	def read(self, path='.', pattern='cyl2flx_diag_*.dat', s_max=None):
		"""Reads and concatenates all per-rank diagnostic files

		Parameters
		----------
		path : str
			Directory containing the diagnostic files.
		pattern : str
			Glob pattern matching the per-rank files.
		s_max : float or None
			If set, keep only points whose solved s (column 4) is below
			this. Use it to drop the deep-vacuum points that dominate a
			run logged with cyl2flx_log_outside=.TRUE. (e.g. s_max=1.5
			focuses on near-plasma / interior failures). None keeps all.

		Returns
		-------
		data : ndarray
			(N, 11) array of all logged points (see module header for
			the column layout). Also stored on self.data.
		"""
		import os
		import glob
		import numpy as np
		self.files = sorted(glob.glob(os.path.join(path, pattern)))
		if not self.files:
			raise FileNotFoundError(
				'No files matching %s in %s. Was LIBSTELL built with '
				'-DDEBUG_CYL2FLX, and did the run produce failures?'
				% (pattern, path))
		blocks = []
		n_raw = 0
		for f in self.files:
			# comments='#' skips the header line; whitespace delimited
			try:
				m = np.genfromtxt(f, comments='#')
			except Exception:
				m = np.empty((0, 11))
			if m.size == 0:
				continue
			m = np.atleast_2d(m)
			if m.shape[1] != 11:
				print('  WARNING: %s has %d columns (expected 11), '
				      'skipping.' % (f, m.shape[1]))
				continue
			n_raw += m.shape[0]
			if s_max is not None:
				m = m[m[:, IS] < s_max]
			blocks.append(m)
		if not blocks or sum(b.shape[0] for b in blocks) == 0:
			raise ValueError('No points left after reading'
			                 + ('' if s_max is None else
			                    ' / s_max filter'))
		self.data = np.vstack(blocks)
		if s_max is None:
			print('Read %d logged points from %d file(s).'
			      % (self.data.shape[0], len(self.files)))
		else:
			print('Read %d points from %d file(s); kept %d with s<%g.'
			      % (n_raw, len(self.files), self.data.shape[0], s_max))
		return self.data

	def classify(self, niter=None, dist_accept=1.0E-3, tau_tol=None,
	             starve_frac=1.0):
		"""Classifies each logged point by likely failure mode

		The classification is ordered by cause so each point gets a
		single, most-informative label:

		  outside    : info == -3  (solver decided s>2; expected for
		               true-exterior points, only present if logged)
		  rescued    : matched within dist_accept -- effectively
		               converged; the caller's fmin_acceptable test
		               accepts these and the field is valid. Tested
		               FIRST because nfe saturates (see note) so the
		               quality of the match, not the effort, is what
		               separates a good point from a bad one.
		  degenerate : NOT matched, and the Jacobian went (near)
		               singular (ndegen>0 or |tau|<=tau_tol)
		                                  -> equilibrium ill-conditioning
		  starved    : NOT matched, ran the full budget without
		               degeneracy            -> tolerance/budget problem
		  failed     : none of the above (genuine non-convergence)

		NOTE on nfe: the value logged is CUMULATIVE over the restart
		loop in cyl2flx (= niter-per-try * number of restart tries, e.g.
		50*8=400), not per-try. Because that loop only exits early on a
		clean hit (fmin<=ftol) or a clear-outside verdict, essentially
		every info=-1 point runs the full budget. nfe therefore saturates
		and is a weak discriminator -- the degeneracy map (taumin/ndegen)
		and the matching error (dist) are the informative quantities.

		Parameters
		----------
		niter : int or None
			Total evaluation budget (cumulative, = per-try cap * number
			of restarts). If None (default) it is inferred from the
			largest nfe present in the dataset, so it tracks the solver
			automatically even if the cap or restart count changes.
		dist_accept : float
			Geometric error [m] below which a point is considered an
			acceptable near-miss (rescued).
		tau_tol : float or None
			|Jacobian| threshold for the degenerate label. If None, the
			same scale newt2d uses is approximated as
			sqrt(eps)*median(R)^2.
		starve_frac : float
			Fraction of the budget at/above which a point counts as
			starved (default 1.0 = used the full budget). Lower it
			(e.g. 0.95) to also catch near-exhausted points.

		Returns
		-------
		mode : ndarray of str
			Per-point mode label (length N).
		"""
		import numpy as np
		if self.data is None:
			raise RuntimeError('Call read() before classify().')
		d = self.data
		n = d.shape[0]
		# Infer the cumulative budget from the data unless overridden.
		if niter is None:
			self.nfe_budget = int(np.max(d[:, INFE]))
		else:
			self.nfe_budget = int(niter)
		if tau_tol is None:
			eps0 = np.sqrt(np.finfo(float).eps)
			tau_tol = eps0 * np.median(d[:, IR])**2
		mode = np.empty(n, dtype=object)
		is_outside = d[:, IINFO] <= -2.5
		is_degen = (d[:, INDEGEN] > 0) | (d[:, ITAUMIN] <= tau_tol)
		is_starved = d[:, INFE] >= starve_frac * self.nfe_budget
		is_rescued = d[:, IDIST] <= dist_accept
		for i in range(n):
			if is_outside[i]:
				mode[i] = 'outside'
			elif is_rescued[i]:
				mode[i] = 'rescued'
			elif is_degen[i]:
				mode[i] = 'degenerate'
			elif is_starved[i]:
				mode[i] = 'starved'
			else:
				mode[i] = 'failed'
		self.mode = mode
		self.tau_tol = tau_tol
		return mode

	def summary(self, niter=None, dist_accept=1.0E-3, tau_tol=None,
	            starve_frac=1.0):
		"""Prints a text summary of failure modes and key statistics

		Parameters
		----------
		niter, dist_accept, tau_tol, starve_frac :
			Passed through to classify(). niter=None infers the budget
			from the data.
		"""
		import numpy as np
		mode = self.classify(niter=niter, dist_accept=dist_accept,
		                      tau_tol=tau_tol, starve_frac=starve_frac)
		d = self.data
		n = d.shape[0]
		print('')
		print('================ cyl2flx convergence summary ============')
		print(' total logged points : %d' % n)
		print(' tau degeneracy tol   : %.3e  (|Ru*Zs-Zu*Rs| <= this)'
		      % self.tau_tol)
		print(' dist_accept          : %.3e m' % dist_accept)
		print(' nfe budget (inferred): %d  (cumulative = niter*restarts)'
		      % self.nfe_budget)
		print(' --- failure modes ---')
		for lab in MODE_LABELS:
			c = int(np.count_nonzero(mode == lab))
			if c > 0:
				print('   %-11s : %7d  (%5.1f%%)'
				      % (lab, c, 100.0 * c / n))
		print(' --- residual / budget ---')
		print('   nfe   median/max  : %d / %d'
		      % (int(np.median(d[:, INFE])), int(d[:, INFE].max())))
		print('   frac at budget    : %5.1f%%'
		      % (100.0 * np.count_nonzero(
		         d[:, INFE] >= self.nfe_budget) / n))
		print('   dist  median/max  : %.3e / %.3e m'
		      % (np.median(d[:, IDIST]), d[:, IDIST].max()))
		print('   frac ndegen>0     : %5.1f%%'
		      % (100.0 * np.count_nonzero(d[:, INDEGEN] > 0) / n))
		print('=========================================================')
		# A one-line verdict to orient the user
		nd = np.count_nonzero((mode == 'degenerate'))
		nstv = np.count_nonzero((mode == 'starved'))
		if nd > nstv and nd > 0:
			print(' VERDICT: degeneracy-dominated -> equilibrium is '
			      'ill-conditioned where it fails (look at the taumin '
			      'map). Tightening tolerance will not help.')
		elif nstv > 0:
			print(' VERDICT: non-degenerate far-misses dominate. Since '
			      'nfe saturates at the budget for every failure, check '
			      'the dist map: if errors are small -> relax the '
			      'accept tolerance; if large -> raise niter/relax '
			      'damping in newt2d, or warm-start.')
		print('')
		return mode

	def plot(self, niter=None, dist_accept=1.0E-3, tau_tol=None,
	         starve_frac=1.0, phi=None, phi_tol=1.0E-3, boundary=None,
	         wout=None, phi_deg=None, show=True):
		"""Produces a 6-panel diagnostic figure

		Parameters
		----------
		niter, dist_accept, tau_tol :
			Passed through to classify().
		phi : float or None
			If given, only points whose Nphi (=nfp*phi) is within
			phi_tol of this value are plotted (select one toroidal
			plane). None plots all planes overlaid.
		phi_tol : float
			Tolerance for the phi selection [rad].
		boundary : tuple(ndarray, ndarray) or None
			Optional (Rb, Zb) of the LCFS to overplot in the (R,Z)
			panels for context.
		wout : str or None
			Path to a VMEC wout file. If given together with phi_deg,
			the LCFS boundary is computed at that angle and overplotted,
			and the matching toroidal plane is auto-selected. This is the
			decisive check for tip failures: it shows whether the failing
			points sit INSIDE or OUTSIDE the truncated-Fourier boundary.
		phi_deg : float or None
			Geometric toroidal angle [degrees] for the wout boundary and
			plane selection.
		show : bool
			Call plt.show() before returning.

		Returns
		-------
		fig : matplotlib Figure
		"""
		import numpy as np
		import matplotlib.pyplot as plt
		if self.data is None:
			raise RuntimeError('Call read() before plot().')
		mode = self.classify(niter=niter, dist_accept=dist_accept,
		                      tau_tol=tau_tol, starve_frac=starve_frac)
		d = self.data
		# If a wout file + angle are given, compute the LCFS boundary at
		# that angle and auto-select the matching plane (snap to nearest
		# logged Nphi so the match is robust to grid spacing).
		if wout is not None and phi_deg is not None:
			Rb, Zb, nfp = lcfs_from_wout(wout, np.deg2rad(phi_deg))
			boundary = (Rb, Zb)
			target = nfp * np.deg2rad(phi_deg)
			phi = float(d[np.argmin(np.abs(d[:, INPHI] - target)),
			              INPHI])
			phi_tol = 1.0E-4
		if phi is not None:
			sel = np.abs(d[:, INPHI] - phi) <= phi_tol
			d = d[sel]
			mode = mode[sel]
			ttl = r' ($N\phi=%.3f$)' % phi
			if d.shape[0] == 0:
				raise ValueError('No points near Nphi=%g.' % phi)
		else:
			ttl = ' (all planes)'

		R, Z = d[:, IR], d[:, IZ]
		fig, ax = plt.subplots(2, 3, figsize=(15, 9))
		fig.suptitle('cyl2flx / newt2d convergence diagnostics' + ttl)

		def _overlay_boundary(a):
			if boundary is not None:
				a.plot(boundary[0], boundary[1], 'k-', lw=1.0,
				       label='LCFS')

		# (1) (R,Z) colored by failure mode
		a = ax[0, 0]
		for lab in MODE_LABELS:
			m = (mode == lab)
			if np.any(m):
				a.scatter(R[m], Z[m], s=10, c=MODE_COLORS[lab],
				          label=lab, edgecolors='none')
		_overlay_boundary(a)
		a.set_title('failure mode')
		a.set_xlabel('R [m]'); a.set_ylabel('Z [m]')
		a.set_aspect('equal', 'box'); a.legend(markerscale=2, fontsize=8)

		# (2) (R,Z) colored by log10(dist)
		a = ax[0, 1]
		sc = a.scatter(R, Z, s=10, c=np.log10(np.maximum(d[:, IDIST],
		               1e-30)), cmap='viridis')
		_overlay_boundary(a)
		fig.colorbar(sc, ax=a, label='log10 dist [m]')
		a.set_title('matching error'); a.set_xlabel('R [m]')
		a.set_ylabel('Z [m]'); a.set_aspect('equal', 'box')

		# (3) (R,Z) colored by log10(taumin) -> degeneracy map
		a = ax[0, 2]
		sc = a.scatter(R, Z, s=10, c=np.log10(np.maximum(d[:, ITAUMIN],
		               1e-30)), cmap='magma')
		_overlay_boundary(a)
		fig.colorbar(sc, ax=a, label='log10 |tau|min [m^2]')
		a.set_title('Jacobian degeneracy'); a.set_xlabel('R [m]')
		a.set_ylabel('Z [m]'); a.set_aspect('equal', 'box')

		# (4) histogram of nfe with niter marker
		a = ax[1, 0]
		a.hist(d[:, INFE], bins=min(50, int(d[:, INFE].max()) + 1),
		       color='tab:blue')
		a.axvline(self.nfe_budget, color='r', ls='--',
		          label='budget=%d' % self.nfe_budget)
		a.set_title('iteration budget'); a.set_xlabel('nfe (cumulative)')
		a.set_ylabel('count'); a.legend(fontsize=8)

		# (5) dist vs Z -> reveals constant-Z banding
		a = ax[1, 1]
		a.semilogy(Z, np.maximum(d[:, IDIST], 1e-30), '.', ms=3,
		           color='tab:blue')
		a.axhline(dist_accept, color='g', ls='--',
		          label='dist_accept=%.0e' % dist_accept)
		a.set_title('matching error vs Z'); a.set_xlabel('Z [m]')
		a.set_ylabel('dist [m]'); a.legend(fontsize=8)

		# (6) mode counts bar chart
		a = ax[1, 2]
		labs = [l for l in MODE_LABELS if np.any(mode == l)]
		cnts = [int(np.count_nonzero(mode == l)) for l in labs]
		a.bar(labs, cnts, color=[MODE_COLORS[l] for l in labs])
		a.set_title('mode counts'); a.set_ylabel('count')
		a.tick_params(axis='x', rotation=30)

		fig.tight_layout(rect=[0, 0, 1, 0.97])
		if show:
			plt.show()
		return fig


def lcfs_from_wout(wout_file, phi_rad, ntheta=361, s_index=-1):
	"""Return the (R,Z) boundary curve of a VMEC equilibrium at one phi.

	Evaluates the last stored flux surface (the LCFS, s_index=-1) of a
	VMEC wout file at the geometric toroidal angle phi_rad, using exactly
	the Fourier coefficients (mpol/ntor) that cyl2flx sees. Overplot this
	on the diagnostic (R,Z) panels to tell whether failing tip points sit
	INSIDE or OUTSIDE the truncated-Fourier boundary:

	  outside the curve -> mode-truncation/ripple; no solver fix helps,
	                       use the extrapolation band (or more VMEC modes)
	  inside  the curve -> genuine basin failure in newt2d

	Requires the compiled LIBSTELL shared library (via libstell.vmec).

	Parameters
	----------
	wout_file : str
		Path to the VMEC wout file.
	phi_rad : float
		Geometric toroidal angle [radians].
	ntheta : int
		Number of poloidal points on the curve.
	s_index : int
		Radial index of the surface (-1 = LCFS).

	Returns
	-------
	Rb, Zb : ndarray
		Boundary R, Z at phi_rad (length ntheta).
	nfp : int
		Number of field periods (maps phi -> Nphi = nfp*phi).
	"""
	import numpy as np
	from libstell.vmec import VMEC
	v = VMEC()
	v.read_wout(wout_file)
	theta = np.linspace(0.0, 2.0 * np.pi, ntheta).reshape(-1, 1)
	phi = np.array([[float(phi_rad)]])
	R = v.cfunct(theta, phi, v.rmnc, v.xm, v.xn)
	Z = v.sfunct(theta, phi, v.zmns, v.xm, v.xn)
	if getattr(v, 'iasym', 0) == 1:
		R = R + v.sfunct(theta, phi, v.rmns, v.xm, v.xn)
		Z = Z + v.cfunct(theta, phi, v.zmnc, v.xm, v.xn)
	Rb = np.asarray(R[s_index, :, 0])
	Zb = np.asarray(Z[s_index, :, 0])
	return Rb, Zb, int(v.nfp)


if __name__ == '__main__':
	import sys
	path = sys.argv[1] if len(sys.argv) > 1 else '.'
	diag = Cyl2FlxDiag()
	diag.read(path=path)
	diag.summary()
	diag.plot()
