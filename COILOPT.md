COILOPT++
=========

This code optimizes a finite set of coils given b-normal on a plasma
surface. \<\<toc\>\>

------------------------------------------------------------------------

### Theory

The COILOPT++ utilizes Levenberg-Marquardt, Genetic, Simulated
Annealing, and Particle Swarm type optimizers to minimize the residual
normal field on a plasma surface. It\'s most widely used mode of
operation is one in which the normal field on an equilibrium boundary is
minimized through variation of coil shapes and surfaces. The plasma
normal field is calculated by the BNORM code. COILOPT++ then utilizes a
B-spline representation of coil traces on a winding surface to calculate
the vacuum field. In it\'s most general form the coil geometry, winding
surface, coil currents, and auxiliary vacuum coil currents may be varied
to minimize the normal field on a plasma boundary. The code is currently
interfaced to the VMEC code.

------------------------------------------------------------------------

### Compilation

This is a C++ code which requires the netCDF, the
[GSL](@http://www.gnu.org/software/gsl/) and
[SILO](@https://wci.llnl.gov/simulation/computer-codes/silo) libraries.
Compilation is straightforward and handled by the include makefile.

------------------------------------------------------------------------

### Input Data Format

The code takes as input a text file which defines the location of other
files and how the code will operate.

    #!C++
    cwsfilename[0] = cws  // Coil winding surface file (multiple may be specified)
    modsplinename[0] = fd.spline  // Modular coil spline file (number matches surface)
    Bmatch.woutname = wout.vmec_test // Text VMEC file name
    Bmatch.cdfname = wout_vmec_test.nc // netCDF VMEC file name
    Bmatch.bnormname = bnorm.vmec_test // Normal field file (as calculated by BNORM code)
    Bmatch.bgridname = bgrid.dat // Normal field file (B on grid?)
    Bmatch.beqoutname = b_norm_eq.dat // Equilibrium normal field grid file (produced by COILOPT++)
    Bmatch.bcoiloutprefix = b_norm_final // Prefix for final coil normal field (produced by COILOPT++)
    tvinname = tv.in
    tvoutprefix = tv
    xyzname = coilxyz.dat
    uvroot = uv
    xyzfized[0] = coils.machine_aux // Auxilliary coils file (one for each auxilliary current group)
    Bmatch.normalize = 1; // On by default.
    Bmatch.nu = 16; // Number of pts in u direction at which to compare
                               // B_normal_plasma with B_normal_coils
                               // Note: if this is left at zero, it will be reset to
                               // 4 times the maximum poloidal mode number in the
                               // Fourier expansion of the plasma field.
    Bmatch.nv = 32; // Number of pts in v direction at which to compare
                                // B_normal_plasma with B_normal_coils
                               // Note: if this is left at zero, it will be reset to
                               // 4 times the maximum toroidal mode number in the
                               // Fourier expansion of the plasma field.

    ncwss = 1;      // Just one winding surface by default.

    Ip_targ = 0.0;  // Use plasma equilibrium value by default.
    curscale[0] = 1.0 // Rescale current by specified amount (number matches surface)
    targlenfactor[0] = 0.0; // coil length to target
    lmodcur[0] = 1 // Vary coil current
    lmodgeom = 1 // Vary coil geometry
    lconstrainIp = 1
    curscal[iws] = 1.0;   // Rescale currents by specified amount
    targlenfactor[iws] = 0.0; // Make coils as short as possible
    lmodcur[iws] = 1;
    lmodgeom = 1;
    lconstrainIp[iws] = 0;
    lsamecur[iws] = 0;
    sadcircvurat[iws] = 0.0;
    curvar = 0.04;
    curautoscal = 1; // Rescale currents to match target R BT at plasma surface
    bgtfR0B0 = 0.0;  // Note: I_TF = 5e+06 * bgtfR0B0
    nfixedgeom = 0;  // By default, don't read in fixed unsplined coilsets.
    lccfixed[0] = 0;    // By default, fixed geom. coil current can vary.
    lverbose = 0;


    nseg_curv = 400;
    nseg_coil = 100;
    dphiguess = 0.1;
    maxcurv = -0.105; // Set negative to use constant coil spacing (faster)
    smartcc = 0;
    minccsep = 0.0; // No penalty for touching coils
    powccsep = 2.0;
    ssexcfracw = 0.1;
      for (int icoil=0; icoil<MAXSC; icoil++) { // No bounds on saddle pos. by default
      sadmodl[icoil] = settings.sadmodr[icoil] = -1;
      sadrectumin[icoil] = settings.sadrectvmin[icoil] = 0.0;
      sadrectumax[icoil] = settings.sadrectvmax[icoil] = 1.0;
      }
    sadmodws = 0;
    sadmodwf = 0.04;
    sadrectpfw = 0.02;

      // chi_sq weights for cost function
    weights.bnorm = 0.0;
    weights.bmax = 1.0;
    weights.Ipol = 0.0;
    weights.curvature = 0.0;
    weights.length = 0.0;
    weights.torsion = 0.0;
    weights.modtorvar = 0.0;
    weights.coilcoil = 0.0;
    weights.selfint = 0.0;
    weights.ssexclusion = 0.0;
    weights.sadmod = 0.0;
    weights.sadrect = 0.0;
    weights.sadcirc = 0.0;

      // DE Optimizer settings (see comments in diffEvol.C for further explanation)
    DE.xtype = rnd; // Use "rnd" to choose random population vectors to breed new
                               // ones to try; use "bst" to use the best current population
                               // vector as the source.
    DE.ndifvecs = 1;  // Number of vectors pairs to subtract to produce new guess
    DE.F  = 0.4;  // Amplification factor for differential variation
    DE.CR = 0.95; // Likelihood of using new vs old source in constructing trial
    DE.atol = 1.0e-6; // Absolute convergence tolerance
    DE.A = 1;DE.B = 4; // Population size = A*(search vector length) + B
    DE.NPmax = 40; // Population upper limit. Use -1 for no upper limit on pop. size
    DE.maxgen = 0;
    DE.savefreq = 100;
    DE.ngenconv = settings.DE.savefreq; // Quit if neither the best nor the worst
                                    // population vector has improved in ngenconv generations.

      // SA Optimizer settings
    SA.llog = 1;  // Turn on logging of annealing progress
    SA.noTmin = 1;  // Turn off minimum temperature convergence criterion
    SA.nvfrac = 0.5;  // Size of subset relative to total vector
    SA.n_uphill_ini = 50; // Number of uphill steps to use in estimating T_0
    SA.rostep = 0.25; // Initial step size as a fraction of total range
    SA.probok = 0.5; // Initial probability of accepting an uphill move
    SA.N1 = 12;  // num_accept < N1*len in temp stage end test
    SA.N2 = 100; // num_try < N2*len in temp stage end test
    SA.rmxtmp = 0.9; // Maximum temperature reduction factor between stages
    SA.rmitmp = 0.1; // Minimum temperature reduction factor between stages
    SA.ratmax = 0.2; // Threshold acceptance ratio for expanding step size
    SA.extstp = 2.0; // Step size expansion factor
    SA.ratmin = 0.05; // Threshold acceptance ratio for shrinking step size
    SA.shrstp = 0.5; // Step size shrinking factor
    SA.nodocr = 0;   // Stopping crit. for successive temp stages w/o downhill mv
    SA.nfmax = 5000; // Maximum fcn evals per vector element
    SA.epsrel = 1.0e-6; // smallest allowed shrink fraction for step vector
    SA.epsabs = 1.0e-8; // smallest allowed value of step vector component

      // PS Optimizer settings
    PS.poprat = 2.75; // Ratio of particle number to degrees of freedom
    PS.omega = 0.5;   // Inertia factor
    PS.phip = 1.0;    // Weight for vector toward personal best
    PS.phig = 1.0;    // Weight for vector toward global best
    PS.convtol = 2.0e-7;   // Mean speed rel to min bounds indicating convergence
    PS.maxiter = 0;        // Max number of iterations to try
    PS.updatefreq = 65536; // How often to generate checkpoint files

      // LM Optimizer settings (library defaults)
      // See lmmin.h for explanation of what these settings control.
    LM.ftol = lm_control_double.ftol;
    LM.xtol = lm_control_double.xtol;
    LM.gtol = lm_control_double.gtol;
    LM.epsilon = lm_control_double.epsilon;
    LM.stepbound = lm_control_double.stepbound;
    LM.maxit = lm_control_double.maxcall;
    LM.printflags = lm_control_double.printflags;
    LM.savefreq = lm_control_double.savefreq;

      // mgrid file settings
    mgrid.lmgrid = 0; // Generate mgrid file?
    mgrid.ir = 51;    // Number of major radius zones in mgrid mesh
    mgrid.jz = 51;    // Number of vertical zones in mgrid mesh
    mgrid.kp = 12;    // Number of toroidal planes in mgrid mesh

      // saddle insertion settings
    saddleinsert.laddsaddle = 0;  // Add a pair of saddle coils? 1 for yes
    saddleinsert.surfid = 0;      // Index of winding surface to add them to
    saddleinsert.ncoefs = 32;     // Number of spline coefficients in coil
    saddleinsert.u0 = 0.2;        // Center u value of coil
    saddleinsert.v0 = 0.2;        // Center v value of coil
    saddleinsert.r_u = 0.05;      // u radius of coil
    saddleinsert.r_v = 0.05;      // v radius of coil
    saddleinsert.current = 0.0;   // counterclockwise coil current in amps
    saddleinsert.lmodgeom = 1;    // Allow new saddle geometry to vary

------------------------------------------------------------------------

### Execution

The code is executed by passing the input parameters file to the code as
such:

    >> mpirun -np 64 ~/bin/xcoilopt++ -f coilopt_params

------------------------------------------------------------------------

### Output Data Format

The code outputs various files based on the type of optimization
preformed, while certain specific files are output no matter the
optimization type.

The equilibrium normal field for the equilibrium, initial configuration,
and final configuration will be output in the b\_norm\_eq.dat,
b\_norm\_init.dat, and b\_norm\_final.dat files. The have similar
formats with the first line specifying the number of toroidal and
poloidal data points, and a series of datapoints. The datapoints are
listed as u, v, b\_normal (where u and v run from 0 to 1)

    76    34
    0.000000e+00    0.000000e+00    0.000000e+00
    1.315789e-02    0.000000e+00    5.538695e-03
    2.631579e-02    0.000000e+00    1.121997e-02
    3.947368e-02    0.000000e+00    1.728472e-02
    5.263158e-02    0.000000e+00    2.369976e-02
    6.578947e-02    0.000000e+00    2.987731e-02
    7.894737e-02    0.000000e+00    3.568746e-02
           .              .              .
           .              .              .
           .              .              .

The final coil spline file is named the same as the initial coil spline
file with \'.out\' appened to the name. The first line contains the
number of field periods, number of coils and number of XXX. Then the
coil data is listed. Each coil is specified by a line beginning with the
number of spline knots, DUU, DUL, the current in the coil, and symmetry.
The next line contains the spline knots (plus 4 additional points). The
U and V values at the spline knot locations are then listed line by line
(number of lines equal to number of knots). The file then contains the
information for the next coil till the end is reached.

    3    4    0
    24    0.000000e+00    0.000000e+00    5.029498e+05    0
    0.0000000e+00 0.0000000e+00 0.0000000e+00 0.0000000e+00 4.7619048e-02 9.5238095e-02 1.4285714e-01 1.9047619e-01 2.3809524e-01 2.8571429e-01 3.3333333e-01 3.8095238e-01 4.2857143e-01 4.7619048e-01 5.2380952e-01 5.7142857e-01 6.1904762e-01 6.6666667e-01 7.1428571e-01 7.6190476e-01 8.0952381e-01 8.5714286e-01 9.0476190e-01 9.5238095e-01 1.0000000e+00 1.0000000e+00 1.0000000e+00 1.0000000e+00
    0.000000e+00    6.275407e-02
    1.587302e-02    6.210337e-02
    2.237615e-02    4.478223e-02
    6.971025e-02    4.783675e-02
    9.728082e-02    2.815039e-02

If requested an mgrid file will be generated (mgrid.lmgrid).

------------------------------------------------------------------------

### Visualization

Explain how to visualize the data.

------------------------------------------------------------------------

### Tutorials

Put links to tutorial pages here.