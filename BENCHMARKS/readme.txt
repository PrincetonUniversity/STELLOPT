This directory contains files for running various benchmarks of the STELLOPT code.


input.STELLOPT_BENCH
This is a basic benchmark which runs each of the non-parallel STELLOPT targets.  Single iteration.


input.STELLOPT_BENCH_TXPORT_GENE
For use when performing a linear GENE run to calculate turbulent transport (serial).  Requires the 'parameters' file. Single iteration.


input.STELLOPT_BENCH_TXPORT_GENE_PARA
For use when performing a linear GENE run to calculate turbulent transport (parallel).  Requires the 'parameters' file. Single iteration.


input.STELLOPT_BENCH_COILOPTPP
For use when performing a COILOPT++ optimization inside STELLOPT.  Requires the 'coilopt_params', 'fd.spline', and 'li383ws.nes' files. (parallel)  250 iterations.