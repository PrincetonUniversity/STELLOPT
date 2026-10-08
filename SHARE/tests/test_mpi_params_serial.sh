#!/bin/sh
# Focused independent behavior oracle; no full-suite serial support claim.
set -eu
solver_repo=$(CDPATH= cd -- "$(dirname "$0")/../.." && pwd)
task_scratch=$(mktemp -d /var/tmp/stellopt-serial-module-test.XXXXXXXX)
trap 'rm -rf "$task_scratch"' EXIT HUP INT TERM
cd "$task_scratch"
cpp -traditional "$solver_repo/LIBSTELL/Sources/Modules/mpi_params.f" mpi_params.f
gfortran -c mpi_params.f
gfortran "$solver_repo/SHARE/tests/test_mpi_params_serial.f90" mpi_params.o -o test_serial
./test_serial > output.log
grep -q 'MPI_STEL_ABORT CALLED BUT NO MPI' output.log
grep -q 'PASS: serial interval and non-MPI diagnostic' output.log
cat output.log
