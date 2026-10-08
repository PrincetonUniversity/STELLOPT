# Native lambda checkpoint normalization

- `wrout.f` exports physical lambda as internal lambda times `lamscale/phipf`, then maps full-grid coefficients to the half grid.
- `load_xc_from_wout.f` reverses the half-grid map and restores `phipf`; it must also restore the inverse `lamscale`.
- `bcovar.f` multiplies internal lambda derivatives by `lamscale`. Neither `initialize_radial.f` nor `profil3d.f` supplies the missing inverse.
- The loader rejects zero, negative or nonfinite normalization with a diagnostic and nonzero exit before reading or division. Finite and positive checks are sequential, so NaN is never compared under floating-point traps.
- The repair changes only restart lambda normalization. Geometry, pressure, iota, interpolation and convergence criteria are unchanged.

Run from the repository root (Python 3 and GNU Fortran):

```sh
python3 tests/test_lambda_restart/run.py --baseline
python3 tests/test_lambda_restart/run.py --require-fixed
```

- The baseline command compiles the loader from official commit `8060f5e5b1bfe11b2f8809c8dfa90459ee72f9ba`; that commit must be present locally.
- Both commands extract the normalization and half-grid export blocks from the actual production `wrout.f` and compile the entire native loader with bounds checking and invalid/division/overflow traps.
- Synthetic read-module data describe exact circular or elongated second-harmonic geometry and physical lambda `0.1 sqrt(s) sin(theta) + 0.2 s sin(2 theta)` on five surfaces.
- Twelve circular/shaped controls use `lamscale = 1, 2, 0.5` and both toroidal-flux signs. `phipf = +/-lamscale` satisfies the native constant-flux normalization; nonzero modes use native `mscale = sqrt(2)`.
- Independent oracles check internal lambda, physical toroidal field density `phipf (1 + d lambda/d theta)`, geometry and pressure/iota invariance.
- Five additional fixed-loader controls must reject zero, negative, NaN and either infinite scale.
- Baseline unit-scale controls pass; its eight nonunit controls must fail the lambda oracle. The fixed loader must pass all twelve with error below `1e-13`.

Scope: axisymmetric stellarator-symmetric component tests, including odd/even poloidal modes and both flux signs. These are not full NetCDF I/O, 3D, asymmetric or multi-rank equilibrium tests. A separate native Make release build checks production integration. No claim is made that this repair resolves cold or warm-start convergence failures.

Startup contract:

- `initialize_radial.f` calls `profil1d` or `profil1d_par` before the only loader call. Both compute `lamscale = sqrt(hs sum(phips(2:ns)^2))` after filling all half-grid flux derivatives.
- With finite nonzero edge toroidal flux and the default linear toroidal-flux law, `lamscale = abs(phiedge)/(2 pi) > 0`. Both flux signs are valid.
- The same-grid `ns_old == ns` early return occurs before either profiles or loader; it cannot expose uninitialized scale at this loader call.
- Existing V3FIT continuation broadcasts its retained profile scale; it does not substitute an inverse normalization.
- The guard also protects direct loader callers. No generic flux-input or initialization redesign is included. Native full-startup checkpoint fidelity remains a separate gate.
