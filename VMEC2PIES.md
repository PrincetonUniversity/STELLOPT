VMEC2PIES
=========

The VMEC2PIES code creates a PIES input file from VMEC input.

------------------------------------------------------------------------

### Theory

The VMEC2PIES code allows a [PIES](PIES) run to be conducted from an
[VMEC](VMEC) input or wout file. This is achieved by calculating a set
of nested surfaces for [PIES](PIES) to use as a background coordinate
system and placing a magnetic field on that grid. It\'s most basic
output is a pies input file with INPUT, PLTFLG, and EXLSTA namelists. It
also outputs the R, Z, B\^S, B\^U, B\^Z, and mode selection matrix in
the input file. The pressure and current profiles are output as a series
of splines over the [VMEC](VMEC) flux coordinates. The code can be
initialize from [VMEC INDATA](VMEC Input Namelist (v8.47)) namelist or a
[VMEC](VMEC) wout file. It can \'extend\' the [VMEC](VMEC) boundary by
extrapolating new surfaces in real space. It can also create a
[PIES](PIES) coil\_data file from a coils file and the EXTCUR array
(found in [INDATA](VMEC Input Namelist (v8.47)) namelist or wout file).
The fixed or free boundary flag in the
[INDATA](VMEC Input Namelist (v8.47)) namelist or wout file can be
overridden for the pies input.

If the user supplies a [VMEC](VMEC) input file with
[INDATA](VMEC Input Namelist (v8.47)) namelist coordinates will be
constructed from the RBC/ZBS arrays and the magnetic axis location
specification. Pressure and current profiles are extracted from the
input file as well. The magnetic field is then calculated by the
following relation: [math](math) B\^u = \\chi\'/\\sqrt(g)\
B\^v = \\Phi\'/\\sqrt(g)\
\\Phi\' = d\\Phi/ds = 1\*\\Phi\\left(s=1\\right)/2\\pi\
\\chi\' = d\\chi/ds = \\iota(s) \* \\Phi\'(s) [math](math) This is the
same initial guess VMEC utilizes for it\'s fields.

If the user supplies a [VMEC](VMEC) wout file, coordinates are generated
by splining over the [VMEC](VMEC) surfaces in real space. Fields are
calculated in a similar manner. External surfaces (if requested) are
extrapolated in real space. The fields in this external region are
either calculated or simply assumed to be the same as found on the
[VMEC](VMEC) boundary. For a free boundary [VMEC](VMEC) run, the vacuum
field is calculated from the mgrid file refrenced in the wout file. The
plasma response is then calculated by a virtual casing principle.

The pressure and current profile (in both cases) will be represented as
a set of spline coefficients over the [VMEC](VMEC) coordinate grid
(normalized toroidal flux). If the total enclosed toroidal flux is zero
the initial scaling factor for [PIES](PIES) (betai) will be set to zero
to guarantee the total current is zero (IOTE will be chosen to be 0 and
ADJST to 0). Otherwise BETAI and IOTE are chosen to match CURTOR for the
current profile supplied.

A coils file (see [MAKEGRID](MAKEGRID)) may also be supplied and the
EXTCUR array from the VMEC file (input or wout) will be utilized to
generate the coil\_data file for [PIES](PIES).

------------------------------------------------------------------------

### Compilation

The VMEC2PIES code has been incorporated in the
[VMEC/STELLOPT](STELLOPT Compilation) makefile system. The code makes
use of the EZSpline libraries
([\@http://w3.pppl.gov/ntcc/PSPLINE/](@http://w3.pppl.gov/ntcc/PSPLINE/))
and LIBSTEL libraries.

------------------------------------------------------------------------

### Input Data Format

The VMEC2PIES code is run from the command line taking the [VMEC](VMEC)
input file or wout file as it\'s first argument

    XVMEC2PIES <VMEC FILE> -n <NTOR> -m <MPOL> -k <SURFS> -extsurfs <EXTERNAL SURFACES> -c <COILS FILE> -fixed -free -noverb -help

\|\| Argument \|\| Default \|\| Description \|\| \|\| -n \|\| VMEC\_FILE
\|\| Maximum selected toroidal mode number \|\| \|\| -m \|\| VMEC\_FILE
\|\| Maximum selected poloidal mode number \|\| \|\| -k \|\| VMEC\_FILE
\|\| Maximum number of surfaces \|\| \|\| -extsurfs \|\| 0 \|\| Number
of PIES surfaces outside the VMEC domain \|\| \|\| -c \|\| NONE\</span\>
\|\| Makegrid style coils file (used to produce PIES coil\_data file, in
conjunction with currents in VMEC File) \|\| \|\| -fixed \|\| \|\|
Overrides the free boundary VMEC flag and forces a fixed boundary PIES
run \|\| \|\| -free \|\| \|\| Overrides the fixed boundary VMEC flag and
forces a free boundary PIES run \|\| \|\| -noverb \|\| \|\| Suppresses
screen output \|\| In it\'s simplest invokation the code requires a VMEC
input or wout file.

------------------------------------------------------------------------

### Execution

The VMEC2PIES code is executed from the command line

    > ~/bin/xvmec2pies wout_test.nc -k 99
     -----PIES File Parameters-----
       extsurfs:  0
              k:  99
     -----VMEC File Parameters-----
        file: test
           m:  11   nu:  45
           n:   0   nv:   1
       mnmax:   12  nuv:   45   nuvp:    45
         nfp:   1
          ns: 385
      lfreeb:   F
        iota: [ 0.107, 1.346]
    torflux_edge: -3.661
    Total Current: 1401793.713
     -----PIES File Parameters-----
        file: test.in
           k:  99 lpinch:  99
           m:  22    mda:  22
           n:   0    nda:   2
         nfp:   1
       freeb:   0
        rmaj:  1.264
       betai: 0.11763476915967E+01
        iote: -.28035874579610E+00
       freeb:   0

------------------------------------------------------------------------

### Output Data Format

The code produces a PIES input file with the INDATA, PLTFLG, EXLSTA
namelists, background coordinates, magnetic fields, mode selection
matrix, pressure, and current splines (see [PIES](PIES)).

------------------------------------------------------------------------

### Visualization

Interested users should see the INIT.rat subroutine in PIES for
interpretation of the quantities in the PIES input file.

------------------------------------------------------------------------

### Tutorials

[NCSX Example](VMEC2PIES NCSX Tutorial)
