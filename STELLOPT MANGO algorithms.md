{% include head.html %} 

The [MANGO library](https://hiddensymmetries.github.io/mango/index.html)
provides STELLOPT with many optimization algorithms from the packages
PETSc, GSL, NLOpt, and HOPSPACK. Many of these algorithms are serial by nature,
and MANGO provides them with parallelization through parallelized finite difference gradient
calculations. Since MANGO is also available from the stellarator optimization system ROSE,
STELLOPT and ROSE can be compared while using exactly the same implementation of the
same optimization algorithm.

An example input file for running STELLOPT with a MANGO algorithm
[can be found here](https://github.com/PrincetonUniversity/STELLOPT/blob/develop/BENCHMARKS/STELLOPT_TEST/MANGO/input.MANGO),
also in `BENCHMARKS/STELLOPT_TEST/MANGO`.


## Building STELLOPT with MANGO
---------------

To use the MANGO algorithms in STELLOPT, you must first build the MANGO library,
following the [directions here](https://hiddensymmetries.github.io/mango/gettingStarted.html).

Next, before compiling STELLOPT, use a text editor to
open the `SHARE/make_*.inc` file for your system. Look for the section on MANGO,
set `LMANGO = T`, and set `MANGO_DIR` appropriately.
Now build STELLOPT following the [standard instructions here](STELLOPT Compilation).


## Selecting a MANGO algorithm
------------------------------

To use an optimization algorithm from the MANGO collection,
just set `OPT_TYPE` to any of the available algorithms
[described here](https://hiddensymmetries.github.io/mango/algorithms.html).
Based on preliminary experience, the most effective algorithms
for STELLOPT problems are the ones designed for least-squares problems
rather than for more general single-objective optimization problems.
Good choices are `gsl_lm`, `gsl_dogleg`, `gsl_ddogleg`, `gsl_subspace2d`, 
`petsc_brgn` (which requires PETSc v3.12 or later), `petsc_pounders`,
and `mango_levenberg_marquardt`.

If your choice of `L*_OPT` and `SIGMA_*` in the `OPTIMUM` namelist means
there are more independent variables than terms in the objective function,
certain algorithms will not work.
In particular, STELLOPT's native Levenberg-Marquardt and GSL's `lm`, `dogleg`, `ddogleg`,
and `subspace2d` algorithms will return an error in this case.
However other algorithms such as `mango_levenberg_marquardt` will still work.


## Additional parameters used by MANGO
------------------------------

Besides `OPT_TYPE`, several other variables in STELLOPT's `OPTIMUM` namelist affect MANGO.

For MANGO algorithms that use finite difference gradients,
the value of `EPSFCN` from the `OPTIMUM` namelist is used for the finite difference step size.

Any bound constraints (a.k.a. box constraints) that are set by the `*_MIN` and `*_MAX` variables
in the `OPTIMUM` namelist are passed to MANGO. These values can be honored by some algorithms but not by others.
Some algorithms *require* bound constraints to be set. For the full list of which algorithms
allow or require bound constraints,
see the [table here](https://hiddensymmetries.github.io/mango/algorithms_8dat_source.html).

`LCENTERED_DIFFERENCES` can be set to `T` or `F` to use either centered or 1-sided
finite differences. This variable only has an effect when `OPT_TYPE` is set to a MANGO algorithm,
i.e. it has no effect on STELLOPT's native Levenberg-Marquardt method, which only uses 1-sided differences.
Note that the setting for `LCENTERED_DIFFERENCES` affects the values of `NOPTIMIZERS` you should choose,
as discussed in more detail below.

The variable `AXIS_INIT_OPTION` affects how the magnetic axis in VMEC is initialized.
Available settings are `"previous"`, `"mean"`, `"midpoint"`, `"input"`, and `"vmec"`.
The default value is `"previous"`. For `"previous"`, the initial guess for the next VMEC calculation will be the             
final axis from the previous VMEC calculation on this processor. This is often an excellent    
guess, but it has the downside that it makes the objective function slightly                   
dependent on the history of previous evaluations of the objective function.                    
This can be a problem when evaluating the finite-difference Jacobian, where                    
small differences in the objective function matter a lot. This history-dependence              
can also introduce dependence on the number of MPI processes, which may be                     
undesirable. It is STRONGLY RECOMMENDED that you set `AXIS_INIT_OPTION` to something
other than `"previous"` to avoid these problems. A recommended setting is `"vmec"`,
which invokes VMEC's internal mechanism for initializing the axis shape as a deterministic
function of the boundary shape. See `STELLOPTV2/Sources/General/stellopt_fcn.f90` for details
of the other options.


## Choosing NOPTIMIZERS
-----------------------

The best value for the `NOPTIMIZERS` variable in STELLOPT's `OPTIMUM` namelist depends
on the number of MPI processes for you job, on your choice of optimization algorithm,
and on your choice of `LCENTERED_DIFFERENCES`.
For algorithms like `petsc_pounders` that do not support concurrent evaluations of the objective function,
you should set `NOPTIMIZERS` to 1. For algorithms like `gsl_dogleg`, `petsc_brgn`, `mango_levenberg_marquardt` etc
that use gradient information, the ideal value of `NOPTIMIZERS` is N + 1 for 1-sided differences,
or 2*N+1 for centered differences, where N is the number of parameters.
Such a value enables the entire finite-difference Jacobian and associated objective function to
be evaluated in the same wallclock time as a single function evaluation.
You may wish to round `NOPTIMIZERS` up to a value that is a factor of the number of processors you will have available, if it is not already a factor.
