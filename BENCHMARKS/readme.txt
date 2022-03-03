This directory contains files for running various benchmarks of the STELLOPT code.
For a full list of options please type

make help

The compare_*.py scripts in these directories depend on the shared library libstell.so (not on libstell.a).
To build this shared library, go to the STELLOPT/LIBSTELL directory and run
make shared_release
You may also need to edit the strings s1, s2, and s3 at the top of STELLOPT/pySTEL/libstell/libstell.py
to match your fortran compiler.
