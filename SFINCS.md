SFINCS is a code that solves the drift-kinetic equation to compute many
neoclassical quantities, such as the collisional radial fluxes, parallel
flows, and bootstrap current. The code has a separate GitHub repository
here <https://github.com/landreman/sfincs> and is not included in the
main STELLOPT repository. This page describes one particular application
of SFINCS: coupling to STELLOPT in order to obtain MHD equilibria with a
self-consistent profile of bootstrap current. For other applications,
you can refer directly to the SFINCS GitHub repository and the user
manual there.

------------------------------------------------------------------------

### Requirements

SFINCS requires a real (not complex) version of
[PETSc](https://www.mcs.anl.gov/petsc/). There is basically no way for
an executable to link to both real and complex versions of PETSc, and
since GENE uses a complex version of PETSc, this means SFINCS and GENE
cannot both be linked to STELLOPT at the same time. (Although perhaps
GENE could be built without PETSc support?)

SFINCS requires that the PETSc library be built with either MUMPS or
SuperLU\_DIST, which are libraries for parallel direct solution of
sparse linear systems. MUMPS is preferable to SuperLU\_DIST, as it is
faster and uses less memory.

SFINCS also requires the HDF5 and NetCDF libraries.

------------------------------------------------------------------------

### Installing, compiling, and linking

To get SFINCS, clone the git repository:

    git clone https://github.com/landreman/sfincs.git

Compile SFINCS following the directions in the SFINCS user manual. A
library file libsfincs.a will be generated that we will link into
STELLOPT.

Next, in your STELLOPT repository, you must (at least for now) switch to
the \"sfincs\" branch of the repo. To do this, type
`git checkout sfincs` from anywhere in the STELLOPT repository.

In your stellopt `make_*.inc` file, look for the section\

    #######################################################################
    # SFINCS Options
    #######################################################################

Here, set `LSFINCS = T` Also, set `SFINCS_DIR` and `PETSC_DIR`
appropriately for your system. Now build STELLOPTV2. This should
generate an executable `xstelloptv2` that is linked to SFINCS.

------------------------------------------------------------------------

### \"vboot\" iterations

To obtain a self-consistent profile of bootstrap current, each call to
(PAR)VMEC is replaced by an iteration between VMEC and a bootstrap
current code (SFINCS, or perhaps some other code like BOOTSJ). VMEC
takes a given profile of toroidal current and computes the magnetic
geometry. Then the bootstrap current code takes the given magnetic
geometry and computes an updated current profile. The process is
repeated until convergence. These iterations are coordinated by
STELLOPT, in the file `STELLOPTV2/Sources/General/stellopt_vboot.f90`.
The switch to turn on this iteration in STELLOPTV2 is
`equil_type = 'vboot'` in the `&optimum` namelist, instead of the usual
settings \'vmec2000\' or \'parvmec\'. To select which code is used to
compute the bootstrap current, you should include
`bootcalc_type='sfincs'` or `bootcalc_type='bootsj'` in the `&optimum`
namelist. In this namelist you should also specify
`vboot_tolerance = 1.e-2` or some other value; this variable controls
the number of iterations that will be performed between VMEC and the
bootstrap current code, with a smaller tolerance leading to more
iterations.

For further mathematical details of the vboot iteration, see the note
`doc/computing_vmec_AC_profile_from_a_bootstrap_current_code`

------------------------------------------------------------------------

### SFINCS input parameters

When SFINCS is run as a standalone code, it reads several namelists
(&general, &geometryParameters, etc.) from a file named input.namelist.
However when SFINCS is run via STELLOPT, these namelists should all be
appended to the main STELLOPT `input.*` file. If there is an
input.namelist file present, it will not be read by SFINCS.

For detailed documentation of the SFINCS input namelists and parameters,
see the SFINCS user manual. The following SFINCS input parameters are
all filled in by STELLOPT, so they do not need to be specified in the
input file, and any values that do appear in the input file will be
over-written:\

    geometryScheme
    VMECRadialOption
    inputRadialCoordinate
    inputRadialCoordinateForGradients
    psiN_wish
    equilibriumFile
    nHats
    THats
    dnHatdpsiNs
    dTHatdpsiNs

There are several parameters related to the SFINCS-STELLOPT interaction
that are not used with standalone SFINCS, and which should be specified
in the `&optimum` namelist. First, `sfincs_s` should be set to a real
array with entries in (0,1\], giving the radii (in terms of normalized
toroidal flux s) at which SFINCS will be run. Do not include 0 as a
radius, since the bootstrap current always vanishes on axis. Depending
on the density and temperature profiles, you may or may not wish to
include 1 as a radius; if either the density or temperature vanish at
s=1 then SFINCS will fail there.

Also, `sfincs_min_procs` should be set to an integer giving the minimum
number of processors allocated to SFINCS for each magnetic surface.

------------------------------------------------------------------------

### Parallelization considerations

The SFINCS calculations at different radii are essentially independent
of each other, and so STELLOPT can efficiently parallelize calculations
over radius. At the same time, SFINCS can take advantage of multiple
processors for the calculation at each radius. An important
consideration is that at each radius, SFINCS requires a substantial
amount of memory (typically 30 GB - 1 TB for experimentally relevant
parameters), and so enough processors must be devoted to each radius or
else SFINCS will crash with an out-of-memory error. Due to this issue,
the STELLOPT parameter `sfincs_min_procs` is provided. STELLOPT will
divide up the processors among radii ensuring that each radius has at
least `sfincs_min_procs` processors. If your SFINCS jobs are running out
of memory but you do not wish to request more processors, raise
`sfincs_min_procs`.

As an example, suppose you are running on the IPP computer Draco or on
NERSC Cori Haswell, either of which has nodes with 32 processors/node
and 128 GB/node, so 4 GB/processor. Suppose further that you are running
SFINCS at a resolution that requires 64 GB per solve, so you need
64/4=16 processors for each instance of SFINCS. Setting
`sfincs_min_procs = 16` is then appropriate. Suppose also that
`sfincs_s` contains 8 radii. If you can request 4 nodes (=4\*32=128
processors), then SFINCS can run for all 8 radii in parallel. If you run
on 2 nodes, then SFINCS will run for radii 1, 3, 5, and 7 in parallel,
and afterwards SFINCS will run for the remaining 4 radii. If you run on
1 node, then SFINCS will first run for radii 1 and 5 in parallel, after
which SFINCS will run for the radii 2 and 6, then radii 3 and 7, then
radii 4 and 8.

There is no restriction that the number of processors must be a multiple
of any particular number related to the number of radii in `sfincs_s` or
to `sfincs_min_procs`. Therefore as long as SFINCS does not run out of
memory, any number of processors should be allowed. For optimal load
balancing, you can request a number of processors equal to the number of
radii in `sfincs_s` times the number of processors you want SFINCS to
use at each radius.

------------------------------------------------------------------------

### SFINCS numerical resolution parameters

It is important to tune the SFINCS resolution parameters Ntheta, Nzeta,
Nxi, and Nx based on the magnetic geometry and collisionality. If these
values are too low, the SFINCS calculation will be under-resolved. If
these resolution parameters are too high, the SFINCS calculations will
take unnecessary computational time and memory. It is recommended that
you do convergence testing using standalone SFINCS (not using STELLOPT)
any time you make significant changes to the density, temperature, or
magnetic geometry. (Small modifications to magnetic geometry, as arise
during optimization, will not require changes to SFINCS resolution, but
major changes such as switching from W7-X to NCSX will require
resolution changes.) For detailed instructions on resolution convergence
testing, see the chapter of the SFINCS user manual on \`Numerical
resolution parameters.\'

------------------------------------------------------------------------

### Radial electric field

The radial electric field Er typically has a small effect on the
bootstrap current, on the order of 10% for W7-X and NCSX. Accurate
computation of Er requires running SFINCS several times at each radius
for various values of Er to solve for the ambipolar value. However,
accurate calculation of Er may not be necessary due to its weak effect
on the bootstrap current. The parameter `sfincs_Er_option` in the
STELLOPT `&optimum` namelist controls how Er is determined in the
STELLOPT-SFINCS system. Presently, two options are available. If
`sfincs_Er_option="zero"`, Er is set to zero. (Actually Er is set to a
very small nonzero value, since SFINCS is slightly more robust when Er
is not exactly zero.) Alternatively, if `sfincs_Er_option="estimate"`,
Er is set to an estimate for the ambipolar radial electric field based
on the sqrt(nu)-regime radial ion flux:

d phi / d s = - (1/e) \[(T\_i / n) (d n / d s) + 1.34 (d T\_i / d s)\]

where phi is the electrostatic potential and s is any radial coordinate.
(This formula is obtained from eq (74.a) of Ho & Kulsrud, Phys Fluids
30, 442 (1987).) Full calculation of the ambipolar Er is not available
yet within STELLOPT, but will be available in the future.

------------------------------------------------------------------------

### Examples

Several example input files are provided in subdirectories `vboot_...`
of the BENCHMARKS directory. The first few lines of each input file
provides some information about the number of processors required,
expected output, etc. Each of the subdirectories also contains a .pdf
file showing the expected convergence of the vboot iterations.

------------------------------------------------------------------------

### Plotting output

Once a vboot calculation (using either BOOTSJ or SFINCS) has run, you
can plot several quantities using the `vbootPlot` script in the `pySTEL`
directory. This script will plot, for instance, how the total toroidal
current and the profile of toroidal current vary over the vboot
iterations. Some of the data this script plots is saved in the
`boot_fit.<extension>` file. The syntax for calling the plotting script
is `vboot_plot boot_fit.<extension> wout_<extension>_vboot0000.nc` You
can add `pdf` at the end of this line to save a .pdf file. Some example
plots are provided in the vboot examples in the BENCHMARKS directory.