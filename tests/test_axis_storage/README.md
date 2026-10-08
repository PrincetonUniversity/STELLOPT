# LASYM axis-storage regressions

- From any directory, run `python3 /path/to/STELLOPT/tests/test_axis_storage/run.py`.
- Run `python3 /path/to/STELLOPT/tests/test_axis_storage/serial_run.py` for the serial path.
- Requires Python3 standard library and GNU Fortran; no MPI, NetCDF, BLAS, full VMEC build or external equilibrium input. `--compiler` selects the GNU-compatible compiler command.
- Each script extracts the actual unchanged native `symforce` routine and matching production save/restore blocks from `funct3d`. The original script checks parallel storage; the serial script checks serial storage. Build files live in a temporary directory and are removed afterward.
- Native force parity outputs overwrite the geometry aliases. Independent sine/cosine oracles check those force values; the source repair restores the original R/Z exactly without changing force arrays.
- Negative control: compile the same native routine/fixture with restoration disabled; the independent geometry-preservation check must fail. The serial test also checks that later iterations leave force scratch untouched. Source hashes and results are printed as JSON.
- Repair scope: LASYM first-step high-force axis repair with `lmove_axis=true`. The serial followup mirrors the reviewed parallel fix. Full native serial solver regressions and multirank behavior are separate gates, not established by this unit.
- Both drivers generate circular and elongated shaped R/Z geometry with a second harmonic. No external equilibrium data are required.
