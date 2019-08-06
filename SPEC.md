SPEC
====

[toc](toc)
[image:SPEC.gif width=\"235\" height=\"360\" align=\"left\"](image:SPEC.gif width="235" height="360" align="left")The
Stepped Pressure Equilibrium Code (SPEC) constructs MHD equilibria by
combining elements of ideal MHD and Taylor relaxation theory in a manner
that is consistent with the structure of chaotic magnetic fields. The
author maintains a detailed documentation site which can be found at
(<http://w3.pppl.gov/~shudson/spec.html>). The key feature of the code
is a stepped pressure profile, which is necessitated as the code only
enforces the topological constraint of flux surface existence at a
finite set of interfaces. Between these interfaces the flat pressure
profile allows for the construction of Beltrami fields on a set of
finite elements. The code is parallelized over the annular regions
between interfaces outputs data in the HDF5 data format.

------------------------------------------------------------------------

### Theory\[\[\#Theory\]\]

Taken from (<http://w3.pppl.gov/~shudson/Spec/descrp.pdf>):

[media type=\"custom\" key=\"13833002\" align=\"center\"](media type="custom" key="13833002" align="center")

------------------------------------------------------------------------

### Compilation\[\[\#Compilation\]\]

The SPEC code is maintained under CVS at PPPL. Please contact the
code\'s author fo access.

------------------------------------------------------------------------

### Input Data Format\[\[\#Input\]\]

The input file format is explained below:
[media type=\"custom\" key=\"13833010\" align=\"center\"](media type="custom" key="13833010" align="center")

Here is a sample input namelist for reference. [code](code) &PHYSICSLIST
LGEOMETRY = 6 LTOROIDAL = 0 LFREEBOUND = 0 PHIEDGE = 1.0 CURTOR = 0.0
EXTCUR(1:100) = 0.0 GAMMA = 1.666666667 NFP = 1.0 NVOL = 8 MPOL = 6 NTOR
= 0 NI(1:8) = 2 2 2 2 2 2 2 2 LCONSTRAINT = 2 TFLUX(1:8) = 0.0 1.0
PFLUX(1:8) = 0.0 1.0 HELICITY(1:8) = 0.0 1.0 PSCALE = 0.0 PRES(1:8) =
0.0 1.0 LADIABATIC = 0 ADIAB(1:8) = 0.0 1.0 BETA = 0.0 MU(1:8) = 0.0 1.0
PL(1:8) = 0 0 QL(1:8) = 1 1 PR(1:8) = 0 0 QR(1:8) = 1 1 MUPFTOL =
1.0E-10 MUPFITS = 1 / &NUMERICLIST IEXTRAP = -1 NDISCRETE = 2 LSLABEL =
0 LSINTERP = 5 NOFE = 2 NQUAD = -1 IMPOL = -1 INTOR = -1 LDUMPF = 0
ISWMIN = 1 LSYMAGL = .FALSE. LPERTURB = 0 DPERTRUB = 1.0E-30 /
&LINEARLIST LBELTRAMI = 1 LINITGUES = 1 LDENSE = 1 LPOSDEF = .TRUE.
SPARSEPC = 1 SPARSEITS = 0 SSOROMEGA = 1.0 SPARSETOL = 0.0 LIOTASOLV = 1
/ &NONLINEARLIST LVACUUM = .FALSE. LC05NDF = .FALSE. LC05PDF = .FALSE.
LE04DGF = .FALSE. LE04LYF = .FALSE. XFTOL = 1.0E-09 FORCEERR = 0.0
FACTOR = 1.0E-02 PWIDTH = 4 QWIDTH = 4 IFINITED = 0 LREADGF = .TRUE.
GFSTUFF = 0 VERIFY = -1 MAXSTEP = -1.0 / &DIAGNOSTICLIST ODETOL =
1.0E-07 ABSREQ = 1.0E-08 RELREQ = 1.0E-08 ABSACC = 1.0E-08 DIVERTORR =
-1.0 DIVERTORZ = 1.0 LPOINCARE = 0 NPPTS = 200 NPTRJ = -1 -1 -1 MPQITS =
10 P1 = 0 0 0 0 P2 = 0 0 0 0 Q1 = 0 0 0 0 Q2 = 0 0 0 0 NPQ = 0 0 0 0
MIRRITS = 0 IRRMPOL = 50 IRRNTOR = 25 IRRSVDCUT = 1.0E-12 IRRTOL =
1.0E-08 LWRPJ = .FALSE. NGHD = 0 LHESSIAN = .FALSE. LHEVALUES = .FALSE.
LHEVECTORS = .FALSE. LCURLERR = .FALSE. LTIMING = .FALSE. / &SCREENLIST
WMA00AB = .FALSE. . . . / &END 0 0 1.5 0.0 1.5 0.0 1.5 0.0 1.5 0.0 1.5
0.0 1.5 0.0 1.5 0.0 1.5 0.0 1 0 0.5 0.5 0.45 0.45 0.4 0.4 0.35 0.35 0.3
0.3 0.25 0.25 0.2 0.2 0.1 0.1 [code](code)

------------------------------------------------------------------------

### Execution\[\[\#Execution\]\]

The SPEC code is executed from the command line via the mpirun command.
The SPEC code takes as input the name of the \'.spec\' file for example
if one wishes to run on 8 processors with the test.spec input file the
command would look like: [code format=\"bash\"](code format="bash") \>
mpirun -np 8 /bin/xspec test [code](code)

------------------------------------------------------------------------

### Output Data Format\[\[\#Output\]\]

The SPEC code outputs various files in the current working directory.
The equilibrium data output file is and HDF5 file. The code also outputs
an input file for restart with the \'.end\' extension. There are also
hidden files (begin with periods) for the poincare plots and the
rotational transform.

------------------------------------------------------------------------

### Visualization\[\[\#Visualization\]\]

Explain how to visualize the data.

------------------------------------------------------------------------

### Tutorials\[\[\#Tutorials\]\]

Put links to tutorial pages here.
