TERPSICHORE
===========

The TERPSICHORE code calculates ideal kink stability from VMEC
equilibria.

TERPSCHORE is managed by the Swiss Plasma Center at the Ecole Polytechnique Federale de Lausanne (EPFL).
STELLOPT only has an interface with it.
If you want to obtain the source code of TERPSCHORE, please contact edith.grueter@epfl.ch or wilfred.cooper@epfl.ch.

After downloading the source code, you can direct the env var `TERPSCHORE_PATH` to the TERPSCHORE folder and set `LTERPSCHORE=T` in the makefile.

------------------------------------------------------------------------

### Theory

The TERPSICHORE code utilizes an ideal MHD model to determine the
stability of the equilibria produced by VMEC. The code performs a
transformation to Boozer coordinates internally (does not make use of
BOOZ_XFORM at this time). The code assumes a perturbation to the VMEC
equilibria of the form 
$$ \vec \xi = \vec \xi \left(
{\xi ^s ,\eta ,\mu } \right) = \sqrt g \xi ^s \nabla \theta
\times \nabla \varphi + \eta
\frac{\vec B \times \nabla s}{B^2} + \left[
{\frac{J\left( s \right)}{\Phi '\left( s \right)B^2 }\eta - \mu
} \right]\vec B $$ 
where the first component is normal to
the flux tube. The perturbation is chosen to be divergence free
eliminating the µ component. Stability is noted by an increase in
potential energy for a given perturbation. Thus negative eigenvalues
indicate unstable modes.

Note: Stability for a mode family is defined by WP/WK being positive.

------------------------------------------------------------------------

### Compilation

The TERPSICHORE code is compiled using a makefile. Before compiling a
code the user must choose the following set of variables and set them
in the `tpr_modules_sp.f` file.

- `NI` Number of VMEC surfaces NS-1
- `IVAC` Number of vacuum surfaces. NS/4 is a good choice
- `NVI=NI+IVAC` Total number of radial grid points
- `MLMNV` The number of equilibrium modes.
- `MLMNB` The number of Boozer modes (first table in input file).
- `MLMNS` The number of stability modes (second table in input file).
- `NJ` Number of poloidal gridpoints used by the code. 
- `NK` Number of toroidal gridpoints to use. 
- `NJK=NJ*NK` Total number of flux surface gridpoints
- `MMAXDF=2*MMAX` where MMAX is the maximum ABS(M) value in the talbes
- `NMAXDF=2*NMAX` where NMAX is the maximum ABS(N) value in the tables
- `ND=NI+IVAC`
- `ND1=ND+1`
- `MD=MLMNS`
- `MDY=MLMNS`
- `NA=2*MLMNS*(NI+IVAC)+MLMNS`

The [pySTEL](pySTEL) utility `terpsichore_util.py` can be used for auto-generating
these numbers from a given VMEC equilibrium.

------------------------------------------------------------------------

### Input Data Format

The TERPSICHORE code is controlled by an input file which is passed to
it via unit 15 (STELLOPT requires this file to be named
terpsichore_input_XX):

                   ARIE3n1
    C
    C        MM  NMIN  NMAX   MMS NSMIN NSMAX NPROCS INSOL
             17    -6    +8    30    -7    11     1     0
    C
    C     TABLE OF FOURIER COEFFIENTS FOR BOOZER COORDINATES
    C     EQUILIBRIUM SETTINGS ARE COMPUTED FROM FIT/VMEC
    C
    C M=  0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6  N
          0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -6
          0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -5
          0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -4
          0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -3
          0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -2
          0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -1
          1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1  0
          1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1  1
          1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1  2
          1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1  3
          1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1  4
          1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1  5
          1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1  6
          1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1  7
          1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1  8
    C
          LLAMPR      LVMTPR      LMETPR      LFOUPR
               0           0           0           0
          LLHSPR      LRHSPR      LEIGPR      LEFCPR
               9           9           1           1
          LXYZPR      LIOTPL      LDW2PL      LEFCPL
               0           1           1           1
          LCURRF      LMESHP      LMESHV      LITERS
               1           1           2           1
          LXYZPL      LEFPLS      LEQVPL      LPRESS
               1           1           0           2
    C
    C    PVAC        PARFAC      QONAX        QN         DSVAC       QVAC    NOWALL
      1.2500e+00  0.0000e-00  0.6500e-00  0.0000e-00  1.0000e-00  1.2500e+00     -2
    C
    C    AWALL       EWALL       DWALL       GWALL       DRWAL       DZWAL   NPWALL
      2.8000e+00  1.8000e+00  5.0000e-01  0.0000e-00  0.0000e-00  0.0000e-00      0
    C
    C    RPLMIN       XPLO      DELTAJP      WCT         CURFAC
      1.0000e-05  1.0000e-06  1.0000e-02  0.0000e-00  1.0000e-00
    C
    C     ANISOTROPIC PRESSURE MODEL APPLIED :                    MODELK =      1
    C
    C     NUMBER OF EQUILIBRIUM FIELD PERIODS PER STABILITY PERIOD: NSTA =      1
    C
    C     TABLE OF FOURIER COEFFIENTS FOR STABILITY DISPLACEMENTS
    C
    C M=  0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6  N
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -7
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -6
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -5
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -4
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -3
          0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -2
          0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  0
          1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  1
          1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  2
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  3
          0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  4
          0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  5
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  6
          1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  7
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  8
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  9
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 11
    C
    C   NEV NITMAX         AL0     EPSPAM IGREEN MPINIT
          1   1500  -5.0 E-03  1.0  E-04      0      0
    C

The following table explains each of these variables [Wiki_TERPS.pdf](docs/Wiki_TERPS.pdf).

| Variable Name | Description |
|----------|:-------------:|
| MM | Maximum poloidal mode number (m=0...MM) in Boozer Spectrum |
| NMIN | Minimum toroidal mode number (n) in Boozer Spectrum |
| NMAX | Maximum toroidal mode number (n) in Boozer Spectrum |
| MMS | Maximum poloidal mode number (m=0...MMS) in Displacement Spectrum |
| NSMIN | Minimum toroidal mode number (n) in Displacement Spectrum |
| NSMAX | Maximum toroidal mode number (n) in Displacement Spectrum |
|NPROCS | Number of processors (not used, should be defaulted to 1)|
| INSOL | 0: VMEC Equilibrium, 1: Solov'ev Equilibrium (radius as radial variable), 2: Solov'ev Equilibrium (volume as radial variable) |
| Boozer Table | This table should match the above spectrum definitions. 0: off, 1:on. |
| LLAMPR | Prints flux surface index i, mode pair index l,m,n, and lambda (16 file 0/1) |
| LVMTPR | Prints the VMEC toroidal angle Boozer Fourier Amplitudes on inner 4 and outer 5 suraces and the Boozer Jacobian amplitudes from 2 alternative reconstructions (16 file, 0/1) |
|LMETPR | Prints the Boozer Fourier amplitudes of R, Z, and VMEC toroidal angle (16 file, 0/1) |
| LFOUPR | NOT USED |
|LLHSPR | Prints the submatrix blocks of the LHS stability matrix and the double Fourier flux tube integrals (16 file, 0/9) |
| LRHSPR | Prints the submatrix blocks of the RHS stability matrix (16 file, 0/9) |
| LEIGPR | NOT USED |
| LEFCPR | NOT USED |
| LXYZPR | NOT USED |
| LIOTPL | NOT USED |
|LDW2PL | NOT USED |
| LEFCPL | Write Xsi and Eta vectors (16 file, 0/1) |
| LCURRF | Controls parallel current density, 1: Reconstructs from charge conservation / MHD force balance, 2: Uses VMEC parallel current density, 9: Construction from metric elements. |
| LMESHP | NOT USED |
| LMESHV | Radial mesh accumulation in the vacuum region, 0 : Exponential, 1 : Equidistant, 2 : Quadratic, 3 : Cubic (recommended), 4 : Quartic towards PVI |
| LITERS | NOT USED |
| LXYZPL | NOT USED |
| LEFPLS | NOT USED |
| LEQVPL | NOT USED |
| LPRESS | Default=0. Otherwise 2 uses VMEC pressure gradient. 9 skips Mercier criterion evaluation |
| PVAC | Exponent governing transition away from PVI to conducting wall (>1) |
| PARFAC | Controls period breaking modes. 0 for periodicity breaking modes. For stllarator symmetry breaking modes (mode number n divisible by number of periods), two modes parities exist 0 and 0.5. |
| QONAX | Q on Axis (for Solov'ev equilibrium) |
| QN | Set ot 0 due to VMEC flux zoning, also applies to TERPSICHORE. |
| DSVAC | Value of radial coordinate s at conducting wall |
| QVAC | Exponent governing transition towards the conducing wall from the PVI (>1) |
| NOWALL | -2: Determine normal at each point of the PVI and rescale by AWALL to obtain conducting wall, (recommended)[Turnbull A D, Cooper W A, Lao L L, and Ku L-P 2011 Ideal MHD spectrum calculations for the ARIES-CS configuration *Nucl. Fusion 51* 123011](http://dx.doi.org/10.1088/0029-5515/51/12/123011), -1 : Conducting wall obtained by multiplying (m/=0) Fourier components by AWALL0 : Conducting wall extrapolated from PVI., 1 : Prescribed conducting wall, Drozdov Formula (GWALL, AWALL, EWALL, DWALL, DRWAL, DZWAL, NPWALL) |
| AWALL | Minor radius of conducting wall. |
| EWALL | Elongation of conducting wall |
| DWALL | Quadrangularity of conducting wall. |
| GWALL | Major Radius of conducting wall |
| DRWAL | Horizontal helical modulation of wall. |
| DZWAL | Vertical helical modulation of conducting wall. |
| NPWALL | Number of toroidal field periods of conducting wall (ignored for NOWALL<1) |
| RPLMIN | Minimum absolute value of R, Z to reprint the active Boozer mode table (6 and 16 file, 1E-5) |
| DELTAJP | Resonance de-tuning parameter for magnetic differential equation (recomend 1E-4 to 0.04) |
| WCT | Horizontal modulation of n=1 m=0 component of wall (nowall=-1). |
| CURFAC | Factor to multiply average parallel curren density in noninteracting fast particle stability model (1.0) |
| MODELK | 0: Noninteracting anisotropic fast particle stability model with reduced kinetic energy, 1: Kruskal-Oberman anisotropic energy principle model with reduced kinetic energy (recommended), 2: Noninteracting anisotropic fast particle stability with physical kinetic energy, 3: Kruskal-Oberman anisotropic energy principle model with physical kinetic energy |
| NSTA | Number of equilibrium periods per stability period (usually equal to equilibrium periods) |
| Displacement Table | This table should match the above spectrum definitions. 0: off, 1: on.|
| NEV | Number of eigenvalue compuations (usually 1, when > 1 it resets AL0 to 95% of previous guess) |
| NITMAX | Number of iterations to converge eigenvalue to that closest to AL0. |
|AL0 | Initial guess for eigenvalue |
| EPSPAM | Relative errof ro eigenvalue convergence. |
| IGREEN | Intended for Green's function solution in vacuum (not implemented) |
| MPINT | The stability mode table is shifted in m by MPINIT. The table usually goes from 0 to 55, with MPIINIT=20 it goes from 20 to 75. |

The equilibrium data is supplied to TERPSICHORE via the fort.18 file.
This file can be generated a few different ways. A python utility is
provided in [pySTEL](pySTEL) `terpsichore_util.py` which can be used to
generate both the equilibrium file and the `terpsichore_input_XX` files
to be passed to the code.  It can be invoked by

    > tpersichore_util.py --vmec VMEC_EXT --input

Note that there appears to be some issue with VMEC's run with `LNYQUIST=T`
so please generate your wout files with `LNYQUIST=F`, which is not the
default VMEC2000 behavior. Also note that this utility will autogenerate
the varialbe list which should be copied and pasted into `tpr_modules_ap.f`.
Only one mode family should be uncommented and the code recompiled for
each run with a different mode family.

------------------------------------------------------------------------

### Execution

The TERPSICHORE code requires two files to run.
 * An equilibrium file named `fort.18`
 * An input file passed on the command line via standard in

The call to TERPSICHORE should look like:

    > xtpr < terpsichore_input_00

------------------------------------------------------------------------

### Output Data Format

The data is output into four files by unit number. 
 1. fort.16 is a log file similar to the screen output but expanded in scope.
 2. fort.17 contains the R and Z values of the input equilibrium grid.
 3. fort.19 contains the plasma boundary and wall harmonics.
 4. fort.22 contains some of the runtime variables.
  * Pressure, Iota, WP, AM(modes), ETA(modes), MU(modes)
 5. fort.23 is a binary file continaing the output dataset.

NOTE: That in this context $$\eta$$ defines the binormal direction.
and $$\mu$$ defines the paralle direction.

------------------------------------------------------------------------

### Visualization

The output files are all text save the `fort.23` file. The [pySTEL](pySTEL)
library can be used for visualizing the data. To do this the user must
compile the TERPSICHORE shared library from the TERPSICHORE source
directory. This is accomplished by `make libterpsichore.so`. Once this
is done and `TERPSICHORE_PATH` is defined the library can be used for
reading and plotting the data of the `fort.23` file.

------------------------------------------------------------------------

### Tutorials

[TERPSICHORE NCSX Tutorial](Tutorial TERPSICHORE NCSX) - old

### References

1. [Anderson D V, Cooper W A, Gruber R, Merazzi S and Schwenn U 1990 Methods for the Efficient Calculation of the (Mhd) Magnetohydrodynamic Stability Properties of Magnetically Confined Fusion Plasmas *International Journal of High Performance Computing Applications* **4** 34--47](http://hpc.sagepub.com/content/4/3/34)

2. [Anderson D V, Cooper W A, Gruber R and Merazzi S 1990 TERPSICHORE: a three-dimensional ideal magnetohydrodynamic stability program *Scientific Computing on Supercomputers*](http://link.springer.com/chapter/10.1007/978-1-4613-0659-7_8)

3. [Fu G Y, Cooper W A, Gruber R, Schwenn U and Anderson D V 1992 Fully three-dimensional ideal magnetohydrodynamic stability analysis of low-n modes and Mercier modes in stellarators *Phys. Fluids B* **4** 1401](http://scitation.aip.org/content/aip/journal/pofb/4/6/10.1063/1.860100)

4. [Redi, M. H., A. Diallo, W. A. Cooper, G. Y. Fu, C. Nührenberg, N. Pomphrey, A. H. Reiman, M. C. Zarnstorff, and NCSX Team 2000 Robustness and flexibility in compact quasiaxial stellarators: Global ideal magnetohydrodynamic stability and energetic particle transport *Physics of Plasmas* **7** 2508](http://dx.doi.org/10.1063/1.874090)
