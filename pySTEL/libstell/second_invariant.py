"""Bounce-action helpers for VMEC; public entry points are VMEC methods.

Uses arrays processed by VMEC.read_wout(): full radial mesh and Fourier
phase m*theta+n*phi. J_over_v includes both bounce legs and has units metres.
"""

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import elementary_charge
from scipy.optimize import brentq, minimize_scalar, minimize


class NoBounceWell(ValueError):
    """No unambiguous trapped well was found at the selected seed point."""


class ForbiddenSeed(NoBounceWell):
    """B(seed) exceeds B_mirror."""


class UnclosedWell(NoBounceWell):
    """A turning point was not found within the finite search window."""


class MarginalOrbit(NoBounceWell):
    """The seed or a turning point is marginal; choose a nearby orbit."""


@dataclass(frozen=True)
class BounceIntegral:
    """Full-bounce reduced action, estimated quadrature error, and bounce angles."""

    J_over_v: float                  # metres: 2 integral sqrt(1-B/B_mirror) dl
    quadrature_difference: float     # metres: successive Gauss-rule difference
    phi_left: float                  # unwrapped geometric toroidal angle
    phi_right: float

    def velocity_invariant(self, energy_eV, mass_kg):
        """Return integral v_parallel dl [m^2/s] for kinetic energy E [eV]."""
        if not (np.isfinite(energy_eV) and energy_eV > 0
                and np.isfinite(mass_kg) and mass_kg > 0):
            raise ValueError("energy_eV and mass_kg must be finite and positive")
        return np.sqrt(2 * energy_eV * elementary_charge / mass_kg) * self.J_over_v

    def momentum_invariant(self, energy_eV, mass_kg):
        """Return integral p_parallel dl [kg m^2/s]. No division by 2*pi."""
        return mass_kg * self.velocity_invariant(energy_eV, mass_kg)


class _Fourier:
    """Paired-point Fourier evaluation: phase m*theta + n*phi after pySTEL read."""

    def __init__(self, m, n, cosine, sine):
        active = (cosine != 0) | (sine != 0)
        self.m, self.n = m[active], n[active]
        self.c, self.s = cosine[active], sine[active]

    def evaluate(self, theta, phi, derivative=False):
        theta, phi = np.broadcast_arrays(theta, phi)
        shape = theta.shape
        t, p = theta.ravel(), phi.ravel()
        value = np.zeros(t.size)
        dtheta = np.zeros(t.size)
        for start in range(0, t.size, 512):
            sl = slice(start, start + 512)
            phase = t[sl, None] * self.m + p[sl, None] * self.n
            co, si = np.cos(phase), np.sin(phase)
            value[sl] = np.sum(co * self.c + si * self.s, axis=1)
            if derivative:
                dtheta[sl] = np.sum((-si * self.c + co * self.s) * self.m, axis=1)
        if derivative:
            return value.reshape(shape), dtheta.reshape(shape)
        return value.reshape(shape)


class _Surface:
    """Fourier coefficients interpolated onto one surface in normalized flux."""

    def __init__(self, vmec, s):
        if not np.isfinite(s) or not 0 < s <= 1:
            raise ValueError("Require 0 < s <= 1; the magnetic axis is excluded")
        ns = int(vmec.ns)
        self.nfp = int(vmec.nfp)
        if ns < 3 or self.nfp < 1:
            raise ValueError("VMEC object has invalid ns or nfp")
        self.iota = float(np.interp(s, np.linspace(0, 1, ns),
                                    np.asarray(vmec.iotaf).ravel()))
        if not np.isfinite(self.iota):
            raise ValueError("Nonfinite VMEC rotational transform")
        asym = bool(int(vmec.iasym))
        x = s * (ns - 1)
        j = min(int(x), ns - 2)
        weight = x - j

        def coeff(name, nmodes):
            a = np.asarray(getattr(vmec, name), dtype=float)
            if a.shape != (ns, nmodes) or not np.all(np.isfinite(a)):
                raise ValueError(f"Invalid {name}; expected finite ({ns}, {nmodes})")
            return (1 - weight) * a[j] + weight * a[j + 1]

        def series(cosname, sinname, nyquist=False):
            suffix = "_nyq" if nyquist else ""
            m = np.asarray(getattr(vmec, "xm" + suffix)).ravel()
            n = np.asarray(getattr(vmec, "xn" + suffix)).ravel()
            c = coeff(cosname, m.size) if cosname else np.zeros(m.size)
            sn = coeff(sinname, m.size) if sinname else np.zeros(m.size)
            return _Fourier(m, n, c, sn)

        self.lam = series("lmnc" if asym else None, "lmns")
        self.b = series("bmnc", "bmns" if asym else None, True)
        self.bsup_phi = series("bsupvmnc", "bsupvmns" if asym else None, True)
        self.r = series("rmnc", "rmns" if asym else None)
        self.z = series("zmnc" if asym else None, "zmns")

    def theta(self, alpha, phi):
        # theta_star = theta + lambda = alpha + iota*phi.
        # phi is the geometric angle, not nfp*phi. xn already includes nfp.
        phi = np.asarray(phi, dtype=float)
        target = (alpha + self.iota * phi + np.pi) % (2 * np.pi) - np.pi
        theta = target.copy()
        p = phi % (2 * np.pi / self.nfp)
        for _ in range(80):
            lam, dlam = self.lam.evaluate(theta, p, derivative=True)
            residual = theta + lam - target
            if np.max(np.abs(residual)) < 2e-12:
                if np.any(1 + dlam <= 0):
                    raise RuntimeError("Non-monotone theta -> theta_star map")
                return theta
            if np.any(np.abs(1 + dlam) < 1e-10):
                raise RuntimeError("Singular theta -> theta_star map")
            theta = theta - np.clip(residual / (1 + dlam), -0.5, 0.5)
        raise RuntimeError("VMEC theta inversion did not converge")

    def field(self, alpha, phi, metric=False):
        theta = self.theta(alpha, phi)
        p = np.asarray(phi) % (2 * np.pi / self.nfp)
        b = self.b.evaluate(theta, p)
        if np.any(~np.isfinite(b)) or np.any(b <= 0):
            raise RuntimeError("Nonpositive or nonfinite reconstructed |B|")
        if not metric:
            return b
        bphi = self.bsup_phi.evaluate(theta, p)
        if np.any(~np.isfinite(bphi)) or np.any(np.abs(bphi) < 1e-14):
            raise RuntimeError("Toroidal parametrization failed: B dot grad(phi) = 0")
        if np.any(bphi > 0) and np.any(bphi < 0):
            raise RuntimeError("B dot grad(phi) changes sign within the bounce well")
        return b, b / np.abs(bphi)  # dl / |dphi|, metres/radian


@lru_cache(maxsize=8)
def _gauss(order):
    return leggauss(order)


def second_adiabatic_invariant(vmec, s, alpha, B_mirror, *, phi_seed=0.0,
                              points_per_period=256, max_periods=32,
                              quadrature_order=64, rtol=1e-7, atol=1e-9):
    """Compute the second invariant for the well containing a section point.

    Parameters
    ----------
    vmec : libstell.vmec.VMEC
        Populated by VMEC.read_wout(). A private _Surface may be reused internally.
    s : float
        Normalized toroidal flux, 0 < s <= 1.
    alpha : float
        Field-line label theta_star - iota(s)*phi, radians.
    B_mirror : float
        E/mu in tesla, using mu = m*v_perp**2/(2*B), E = m*v**2/2.
        Keep this fixed across a contour map at fixed E and mu.
    phi_seed : float
        Unwrapped toroidal angle of the seed. The nearest turning points to
        either side select exactly one connected well; other wells are excluded.
    points_per_period : int
        Samples per 2*pi/nfp in the bounce-point search. Increase to resolve
        narrow ripple barriers. This controls topology, not quadrature accuracy.
    max_periods : int
        Maximum search distance in field periods ON EACH SIDE of phi_seed.
        Failure to close a well is not proof that the orbit is passing.
    quadrature_order : int
        Initial Gauss-Legendre order. Doubled up to three times for convergence.
    rtol, atol : float
        Tolerances on successive quadrature rules for J_over_v (atol in metres).

    Returns
    -------
    BounceIntegral
        J_over_v = 2*integral_left^right sqrt(1-B/B_mirror)*dl [m].
        Multiply by v=sqrt(2*E/m) for integral v_parallel dl; multiply by m*v
        for integral p_parallel dl. The factor two includes the return leg.

    Notes
    -----
    Assumes a static magnetic field, zero electrostatic potential, nonrelativistic
    motion, and negligible orbit width. This is a field-line bounce action;
    finite-width guiding-centre trajectories are not integrated here.
    Surface interpolation is linear in s after pySTEL's half-to-full conversion.
    Resolve and converge bounce-point sampling independently of quadrature.
    """
    if not np.isfinite(B_mirror) or B_mirror <= 0:
        raise ValueError("B_mirror must be finite and positive")
    if not np.isfinite(alpha) or not np.isfinite(phi_seed):
        raise ValueError("alpha and phi_seed must be finite")
    for name, value, minimum in [("points_per_period", points_per_period, 16),
                                  ("max_periods", max_periods, 1),
                                  ("quadrature_order", quadrature_order, 16)]:
        if int(value) != value or value < minimum:
            raise ValueError(f"{name} must be an integer >= {minimum}")
    if not (np.isfinite(rtol) and rtol > 0 and np.isfinite(atol) and atol >= 0):
        raise ValueError("Require finite rtol > 0 and atol >= 0")
    surface = vmec if isinstance(vmec, _Surface) else _Surface(vmec, s)
    period = 2 * np.pi / surface.nfp

    def q(phi):
        return 1 - surface.field(alpha, phi) / B_mirror

    qseed = float(q(phi_seed))
    if qseed < -1e-12:
        raise ForbiddenSeed("Seed is inaccessible: B(seed) > B_mirror")
    if qseed <= 1e-12:
        raise MarginalOrbit("Seed lies at or too close to a turning point")

    def turning_point(direction):
        # An extra point beyond each chunk resolves extrema at chunk boundaries.
        # Refine sampled field maxima as well as sign changes: a shallow, narrow
        # forbidden barrier can sit between two accessible samples.
        offsets = np.arange(int(points_per_period) + 2) / points_per_period
        for k in range(int(max_periods)):
            distance = period * (k + offsets)
            phi = phi_seed + direction * distance
            values = q(phi)
            minima = np.flatnonzero((values[1:-1] <= values[:-2])
                                    & (values[1:-1] < values[2:])) + 1
            extra_distance, extra_values = [], []
            for j in minima:
                minimum = minimize_scalar(
                    lambda d: float(q(phi_seed + direction * d)),
                    bounds=(distance[j - 1], distance[j + 1]), method="bounded",
                    options={"xatol": period * 1e-12})
                if not minimum.success:
                    raise RuntimeError("Magnetic-barrier refinement failed")
                extra_distance.append(minimum.x)
                extra_values.append(minimum.fun)
            distances = np.concatenate([distance[:-1], extra_distance])
            vals = np.concatenate([values[:-1], extra_values])
            marginal = np.concatenate([np.zeros(distance.size - 1, dtype=bool),
                                       np.abs(extra_values) <= 1e-12])
            order = np.argsort(distances)
            distances, vals = distances[order], vals[order]
            marginal = marginal[order]
            barrier = np.flatnonzero((vals <= 0) | marginal)
            if barrier.size:
                j = barrier[0]
                if marginal[j]:
                    raise MarginalOrbit("Magnetic barrier is too close to a separatrix")
                if j == 0:
                    raise RuntimeError("Unresolved barrier at search-chunk boundary")
                a = phi_seed + direction * distances[j - 1]
                b = phi_seed + direction * distances[j]
                root = brentq(lambda p: float(q(p)), min(a, b), max(a, b),
                              xtol=1e-12, rtol=1e-13)
                h = period * 1e-5
                slope = abs(float(q(root + h) - q(root - h))) / (2 * h)
                if slope * period < 1e-8:
                    raise MarginalOrbit("Turning point is too close to a separatrix")
                return root
        raise UnclosedWell(f"No turning point in direction {direction:+d} within "
                           f"{max_periods} field periods; extend the search")

    left, right = turning_point(-1), turning_point(+1)
    midpoint, halfwidth = (left + right) / 2, (right - left) / 2

    def integrate(order):
        # phi = midpoint + halfwidth*sin(t) removes the endpoint square roots.
        nodes, weights = _gauss(order)
        t = (np.pi / 2) * nodes
        phi = midpoint + halfwidth * np.sin(t)
        b, dl_dphi = surface.field(alpha, phi, metric=True)
        radicand = 1 - b / B_mirror
        if np.min(radicand) < -1e-10:
            raise RuntimeError("Unresolved internal magnetic barrier: increase "
                               "points_per_period")
        jacobian = halfwidth * np.cos(t) * np.pi / 2
        return float(2 * np.sum(weights * np.sqrt(np.maximum(radicand, 0))
                                * dl_dphi * jacobian))

    value = integrate(int(quadrature_order))
    for multiplier in (2, 4, 8):
        refined = integrate(int(quadrature_order) * multiplier)
        difference = abs(refined - value)
        if difference <= atol + rtol * abs(refined):
            return BounceIntegral(refined, difference, left, right)
        value = refined
    raise RuntimeError("Bounce quadrature did not converge; increase quadrature_order")


def invariant_grid(vmec, s_values, alpha_values, B_mirror, **options):
    """Return (J_over_v, status) arrays shaped (n_s, n_alpha).

    Status: 'ok', 'forbidden', 'unclosed', or 'marginal'. Undefined values are
    NaN, never zero. Runtime/numerical errors propagate instead of being hidden.
    Supply the same phi_seed for the whole section. This labels wells by their
    intersection with that section, not by a sorted well index on a long line.
    """
    s_values = np.asarray(s_values, dtype=float)
    alpha_values = np.asarray(alpha_values, dtype=float)
    if s_values.ndim != 1 or alpha_values.ndim != 1:
        raise ValueError("s_values and alpha_values must be one-dimensional")
    if not s_values.size or not alpha_values.size:
        raise ValueError("s_values and alpha_values must be nonempty")
    result = np.full((s_values.size, alpha_values.size), np.nan)
    status = np.full(result.shape, "unclosed", dtype="U10")
    for i, s in enumerate(s_values):
        surface = _Surface(vmec, s)
        for j, alpha in enumerate(alpha_values):
            try:
                bounce = second_adiabatic_invariant(surface, s, alpha,
                                                   B_mirror, **options)
            except ForbiddenSeed:
                status[i, j] = "forbidden"
            except UnclosedWell:
                status[i, j] = "unclosed"
            except MarginalOrbit:
                status[i, j] = "marginal"
            else:
                result[i, j] = bounce.J_over_v
                status[i, j] = "ok"
    return result, status


@dataclass(frozen=True)
class FieldExtrema:
    """Refined |B| extrema [T] and their (s, theta_VMEC, geometric phi) locations."""

    b_min: float
    b_max: float
    s_min: float
    theta_min: float
    phi_min: float
    s_max: float
    theta_max: float
    phi_max: float


def field_extrema(vmec, ntheta=64, nphi=64):
    """Find equilibrium-wide magnetic mirror-field bounds and the seed section.

    Scan every full-grid surface, including axis and boundary, over one field
    period. Angular sampling is at least four points per highest harmonic;
    refine each surface's sampled minimum and maximum with analytic gradients.
    These are numerical extrema, not rigorous bounds for an unresolved spectrum.
    With linear interpolation of coefficients in s, global extrema occur on
    full-grid surfaces. The axis is safe here because no line metric is used.
    """
    m = np.asarray(vmec.xm_nyq, dtype=float).ravel()
    n = np.asarray(vmec.xn_nyq, dtype=float).ravel() / vmec.nfp
    if int(ntheta) != ntheta or int(nphi) != nphi or min(ntheta, nphi) < 8:
        raise ValueError("ntheta and nphi must be integers >= 8")
    ntheta = max(int(ntheta), int(4 * np.max(np.abs(m))) + 1)
    nphi = max(int(nphi), int(4 * np.max(np.abs(n))) + 1)
    theta = np.linspace(0, 2 * np.pi, ntheta, endpoint=False)
    v = np.linspace(0, 2 * np.pi, nphi, endpoint=False)
    c = np.asarray(vmec.bmnc, dtype=float)
    sn = (np.asarray(vmec.bmns, dtype=float) if vmec.iasym
          else np.zeros_like(c))
    b = vmec.cfunct(theta[:, None], v[:, None] / vmec.nfp,
                    c, vmec.xm_nyq, vmec.xn_nyq)
    if vmec.iasym:
        b += vmec.sfunct(theta[:, None], v[:, None] / vmec.nfp,
                         sn, vmec.xm_nyq, vmec.xn_nyq)
    if not np.all(np.isfinite(b)) or np.min(b) <= 0:
        raise RuntimeError("Nonpositive or nonfinite equilibrium |B|")
    extrema = []
    for sign in (1, -1):
        best = None
        for i, surface in enumerate(b):
            j, k = np.unravel_index(np.argmin(sign * surface), surface.shape)
            initial = np.array([theta[j], v[k]])

            def objective(angles):
                phase = m * angles[0] + n * angles[1]
                co, si = np.cos(phase), np.sin(phase)
                value = np.sum(c[i] * co + sn[i] * si)
                derivative = -c[i] * si + sn[i] * co
                gradient = np.array([np.sum(derivative * m),
                                     np.sum(derivative * n)])
                return sign * value, sign * gradient

            result = minimize(objective, initial, method="BFGS", jac=True,
                              options={"gtol": 1e-10, "maxiter": 100})
            value, _ = objective(initial)
            angles = initial
            # BFGS can report precision loss at an otherwise converged extremum.
            if np.isfinite(result.fun) and result.fun <= value:
                value, angles = float(result.fun), result.x
            if best is None or value < best[0]:
                angles = np.mod(angles, 2 * np.pi)
                angles[np.isclose(angles, 2 * np.pi, rtol=0, atol=1e-10)] = 0
                best = (float(value), i / (vmec.ns - 1),
                        float(angles[0]), float(angles[1] / vmec.nfp))
        extrema.append(best)
    lo, hi = extrema
    return FieldExtrema(lo[0], -hi[0], *lo[1:], *hi[1:])


@dataclass(frozen=True)
class JInvariantGrid:
    """Reduced full-bounce action on a rho,alpha grid; undefined actions are NaN."""

    rho: np.ndarray
    alpha: np.ndarray
    J_over_v: np.ndarray
    status: np.ndarray
    b_mirror: float
    phi_seed: float
    extrema: FieldExtrema


def invariant_rho_grid(vmec, b_mirror=4.9, rho=None, nalpha=128,
                       phi_seed=None, extrema=None, **options):
    """Compute J/v on rho=sqrt(s), alpha=theta_star-iota*phi.

    The default radius grid has 64 points from the first non-axis VMEC radius
    through the LCFS. Auto-seeding uses the geometric phi of the global B minimum
    as one common section, preserving the same well-selection rule across the
    map. Other wells may not intersect that section. Alpha has nalpha distinct
    samples; plotting separately closes the periodic seam.
    """
    if not np.isfinite(b_mirror) or b_mirror <= 0:
        raise ValueError("b_mirror must be finite and positive")
    if int(nalpha) != nalpha or nalpha < 4:
        raise ValueError("nalpha must be an integer >= 4")
    if rho is None:
        rho = np.linspace(np.sqrt(1 / (vmec.ns - 1)), 1, 64)
    rho = np.asarray(rho, dtype=float)
    if (rho.ndim != 1 or rho.size < 2 or not np.all(np.isfinite(rho))
            or np.any(rho < 0) or np.any(rho > 1) or np.any(np.diff(rho) <= 0)):
        raise ValueError("rho must increase strictly within [0,1], with >= 2 points")
    if extrema is None:
        extrema = field_extrema(vmec)
    if phi_seed is None:
        phi_seed = extrema.phi_min
    if not np.isfinite(phi_seed):
        raise ValueError("phi_seed must be finite")
    alpha = np.linspace(0, 2 * np.pi, int(nalpha), endpoint=False)
    j = np.full((rho.size, alpha.size), np.nan)
    status = np.full(j.shape, "axis", dtype="U10")
    positive = rho > 0
    if b_mirror < extrema.b_min:
        status[positive] = "forbidden"
    elif b_mirror >= extrema.b_max:
        status[positive] = "unclosed"
    else:
        j[positive], status[positive] = invariant_grid(
            vmec, rho[positive] ** 2, alpha, b_mirror,
            phi_seed=phi_seed, **options)
    return JInvariantGrid(rho, alpha, j, status, float(b_mirror),
                          float(phi_seed), extrema)


def plot_invariant_rho(grid, ax=None, levels=20, cmap="viridis"):
    """Draw a periodic polar contour map and return (figure, polar axes, grid).

    cmap accepts a Matplotlib colormap name or Colormap (default: viridis).
    """
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 7), subplot_kw={"projection": "polar"},
                                constrained_layout=True)
    else:
        if ax.name != "polar":
            raise ValueError("ax must use a polar projection")
        fig = ax.figure
    alpha = np.append(grid.alpha, grid.alpha[0] + 2 * np.pi)
    values = np.ma.masked_invalid(np.column_stack((grid.J_over_v, grid.J_over_v[:, 0])))
    if np.isfinite(grid.J_over_v).sum() >= 4:
        color = ax.contourf(alpha, grid.rho, values, levels=levels,
                            cmap=cmap, corner_mask=False)
        if np.ptp(values.compressed()) > 0:
            ax.contour(alpha, grid.rho, values, levels=color.levels,
                        colors="k", linewidths=1, alpha=.65, corner_mask=False)
        fig.colorbar(color, ax=ax, pad=.10,
                      label=r"$J_\parallel/v$ [m] (full bounce)")
    else:
        ax.text(.5, .5, "No closed bounce wells on this grid",
                 transform=ax.transAxes, ha="center", va="center")
    ax.set_ylim(0, 1)
    ax.set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2],
                  ["0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$"])
    ax.set_xlabel(r"$\alpha=\theta_* - \iota\varphi$ [rad]")
    ax.set_ylabel(r"$\rho=\sqrt{s}$", labelpad=25)
    ax.set_title("Second adiabatic invariant\n"
                  + rf"$B_{{\rm mirror}}={grid.b_mirror:g}$ T, "
                  + rf"$\varphi_{{\rm seed}}={grid.phi_seed:.6g}$ rad", pad=24)
    return fig, ax, grid
