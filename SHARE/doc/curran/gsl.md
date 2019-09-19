Binding to C optimizers
=======================

Quick start
-----------

1. Check out STELLOPT and set up Docker container
   
   Use `SHARE/docker-build.sh` for container creation and first build, and use `docker-run.sh` for opening a terminal.

2. Build `xstelloptv2`
   
   Invoke `make` from the "STELLOPTV2" directory (or run `./build-all` from the root directory).  (If you just created a fresh container with `docker-build.sh`, this step is unnecessary.)

3. Check out "stellopt\_scenarios" repository
   
   `git clone git@github.com:landreman/stellopt_scenarios.git`.  The steps below assume this is checked out as a sibling to the STELLOPT repository inside the container.

4. Modify input file to request GSL optimizer
   
   Edit 1DOF\_circularCrossSection\_varyAxis\_targetIota/Levenberg-Marquardt/input.test", uncomment line 52, and replace its contents with `opt_type = 'clmdif'`.  `vi` (vim-tiny) is included in the container.

5. Run STELLOPT
   
   From the "1DOF\_circularCrossSection\_varyAxis\_targetIota/Levenberg-Marquardt" directory, run `../../../STELLOPT/STELLOPTV2/Release/xstelloptv2 input.test` (this assumes that the STELLOPT and stellopt\_scenarios directories are siblings).
   
   For each evaluation of the objective function, `ncnt` will be printed.  For each iteration of the minimizer, `iter`, `x`, and `f` will be printed.
   
   Note that `mpirun` should NOT be used when executing `xstelloptv2`, as the C layer is not currently MPI-aware.

TODO
----

* MPI support

GSL wrapper notes
-----------------

* GSL strives to give users lots of control over how numerical algorithms are composed and executed.  However, it doesn't always provide enough control over system-level details like memory allocation.  In particular, the minimizers allocate their own arrays for `x` and `f`, the latter of which may be a slice of a Jacobian matrix with non-unit stride.  This makes it awkward to take advantage of the storage for these variables that has already been allocated on the Fortran side.
  Currently we pass the arrays provided by GSL to the objective function callback when possible and only copy back into the Fortran arrays at the end.  But if the stride of `f` is not 1, then we use the Fortran array `fvec` and copy out of it for each function evaluation.
* GSL's finite difference Jacobian calculation is too simple to be robust with objective functions like STELLOPT's.  Currently we pass through a parameter meant to adjust the step size and set GSL's absolute `h` to the square root of that value.  However, MINPACK actually treats that as a relative step size, so a hack may be necessary to achieve rough equivalence.  For example, scaling `h` by 0.1 is appropriate for the 2DOF example case.
* This wrapper was written at midnight; there may be bugs.

Language binding basics
-----------------------

When invoking code written in one language from another, there are some low-level details that need to be controlled:

* Calling conventions
  
  While the idea of a "function" seems pretty universal across programming languages, the act of calling one relies on conventions that may differ even between compilers for the same language.  These "calling conventions" include which arguments get passed in registers vs. on the stack, where the return value is stored, etc.
  
  To call a C function from Fortran, Fortran needs to know to use C's calling conventions for that invocation.  This is generally controlled by compiler flags; when using the GNU compiler collection with default flags, conventions will typically match.

* Symbol naming
  
  The symbols that object files make available for linking often have names that have been modified or mangled relative to what you type in the source code (to support overloading, for instance, or to canonicalize case-insensitive names).  Different languages have different conventions, so in order to allow them to find each other's symbols at link time, you may need to change the default name generation scheme (or possibly tweak names manually).
  
  In Fortran, `BIND(C)` (optionally with the `NAME` specifier) yields C-compatible names.  In C++, `extern C` accomplishes the same.

* Type representation
  
  Integers may have different sizes, and floating-point numbers may have different precisions.  Struct fields may be ordered and padded differently.  Basic compound data types (strings, multi-D arrays) may have different conventions.
  
  In Fortran 2003, the `ISO_C_BINDING` module provides `KIND` values matching primitive C types on the platform.  For higher-level types, do your homework (or decompose into primitives).

* Argument passing
  
  By default, Fortran passes all arguments by reference, while C typically passes primitive types by value.  This can be controlled in Fortran with attributes like `VALUE` and `ATTRIBUTES C`, or C function signatures can be changed to expect pointers for all parameters.

* Aliasing
  
  Most relevant when calling Fortran from C, but has performance implications the other way around.  See `restrict`.

To pass a callback from one language to another, most of the above (except naming) must be considered for the callback function as well.  A Fortran subroutine can be passed as a C function pointer using `C_FUNLOC`.

References
----------

* [GSL Nonlinear Least-Squares Fitting](https://www.gnu.org/software/gsl/doc/html/nls.html)
* [Introduction to Modern Fortran: Interoperability with C](http://people.ds.cam.ac.uk/nmm1/Fortran/paper_14.pdf)
