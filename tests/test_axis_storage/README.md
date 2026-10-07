# LASYM axis-storage regressions

- From any directory, run `python3 /path/to/STELLOPT/tests/test_axis_storage/run.py`.
- Run `python3 /path/to/STELLOPT/tests/test_axis_storage/serial_run.py` for the serial path.
- Requires Python3 standard library and GNU Fortran; no MPI, NetCDF, BLAS, full VMEC build or external equilibrium input. `--compiler` selects the GNU-compatible compiler command.
- Each script extracts the actual unchanged native `symforce` routine and matching production save/restore blocks from `funct3d`. The original script checks parallel storage; the serial script checks serial storage. Build files live in a temporary directory and are removed afterward.
- Native force parity outputs overwrite the geometry aliases. Independent sine/cosine oracles check those force values; the source repair restores the original R/Z exactly without changing force arrays.
- Negative control: compile the same native routine/fixture with restoration disabled; the independent geometry-preservation check must fail. The serial test also checks that later iterations leave force scratch untouched. Source hashes and results are printed as JSON.
- Repair scope: LASYM first-step high-force axis repair with `lmove_axis=true`. The serial followup mirrors the reviewed parallel fix. Full native serial solver regressions and multirank behavior are separate gates, not established by this unit.
- Additional retained native checks: a symmetric analytical E0 `wout` is byte-identical before/after the patch; an unchanged-input TC24 cold start changes from nonphysical axis/NaN/ier16 to physical axis/native convergence. These larger campaign inputs are separate from this portable unit.
- Native convergence does not establish full physical accuracy. TC24 continuous force and independent resolution refinement remain open.
