BNORM
=====

![](images/bnorm_ncsx_c09r00_free.png)

The BNORM code calculates the Fourier coefficients of the magnetic field
normal to the VMEC surface. [math](math)
\\vec{B}\\cdot\\hat{n}=\\frac{1}{4\\pi}\\int\\frac{\\vec{j}\'
dA}{\|\\vec{x}-\\vec{x}\'\|\^{3/2}} [math](math) It achieves this by
utilizing the covariant magnetic field on the VMEC surface to construct
a surface current. It creates input files for the NESCOIL program.

------------------------------------------------------------------------

### Theory

The purpose of this code is to calculate the Fourier components of
outward normal components of the magnetic field. This information is
passed to the NESCOIL routine. The code achieves this by first using the
covariant magnetic field to construct a surface current density
[math](math) \\vec{K}=B\_u \\frac{\\partial \\vec{x}}{\\partial v} -
B\_v \\frac{\\partial \\vec{x}}{\\partial u}. [math](math) This surface
current is then used to construct the vector potential everywhere on
VMEC surface. The calculation of the vector potential requires special
treatment for singularities. At a given point x the vector potential due
to all other points on the surface is calculated. Then a subroutine is
called to handle contribution to the vector potential by the point x.
The curl of this vector potential then provides the normal component of
the magnetic field on the VMEC surface.

------------------------------------------------------------------------

### Compilation

BNORM is a component of the STELLOPT suite of codes.

------------------------------------------------------------------------

### Input Data Format

The BNORM code takes the VMEC \'wout\' file as it\'s only input file.

------------------------------------------------------------------------

### Execution

To run BNORM with a given input file simply pass the name of the VMEC
output file to BNORM like so (input file named input.test):

    yourmachine:0005> ~/bin/xbnorm wout.test >& bnorm_log.test &
    or
    yourmachine:0005> ~/bin/xbrnom wout.test 0.3 >& bnorm_log.test &

Here we\'ve redirected screen output (trapping error messages) to
\'bnorm\_log.test\' and put the process in the background. The code
takes one optional parameter, the coil separation in generalized units.

------------------------------------------------------------------------

### Output Data Format

The BNORM code outputs two files the \'bnorm.suffix\' file and
\'nescin.suffix\' file (where suffix is the suffix of the VMEC \'wout\'
file). The first file contains the Fourier coefficients of the normal
B-Field. The second file is the input file for the NESCOIL code.

The \'bnorm\' file has a list of the m, n, and bnormal (sine) Fourier
coefficients. For information on the \'nescin\' file please see the
NESCOIL page.

------------------------------------------------------------------------

### Visualization

There are no specific routines to visualize the data in the \'bnorm\'
file. However, the matlabVMEC routines should allow a user to plot the
values stored in the \'bnorm\' file. The bnorm values are calculated
from the covariant components of the magnetic field from VMEC. Thus, the
Fourier coefficients should have the same parity (sine) as the covariant
components.

------------------------------------------------------------------------

### Tutorials

After running VMEC, simply call the code with the VMEC \'wout\' file as
it\'s input parameter as explained above. This should produce the
\'bnorm\' and \'nescin\' files.
