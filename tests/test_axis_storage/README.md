# Parallel LASYM axis-storage regression

- From any directory, run `python3 /path/to/STELLOPT/tests/test_axis_storage/run.py`.
- Requires Python3 standard library and GNU Fortran; no MPI, NetCDF, BLAS, full VMEC build or external equilibrium input. `--compiler` selects the GNU-compatible compiler command.
- The script extracts the actual unchanged `symforce_par` routine and the production save/restore blocks from `funct3d_par`. Build files live in a temporary directory and are removed afterward.
- Native force parity outputs overwrite the geometry aliases. Independent sine/cosine oracles check those force values; the source repair restores the original R/Z exactly without changing force arrays.
- Negative control: compile the same native routine/fixture with restoration disabled; the independent geometry-preservation check must fail. Source hashes and both results are printed as JSON.
- Repair scope: parallel LASYM first-step high-force axis repair. Serial-path repair and multirank behavior are not established by this unit.
- Additional retained native checks: a symmetric analytical E0 `wout` is byte-identical before/after the patch; an unchanged-input TC24 cold start changes from nonphysical axis/NaN/ier16 to physical axis/native convergence. These larger campaign inputs are separate from this portable unit.
- Native convergence does not establish full physical accuracy. TC24 continuous force and independent resolution refinement remain open.
