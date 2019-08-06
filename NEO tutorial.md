Tutorial: Running NEO for NCSX
==============================

In this tutorial we will run NEO as a standalone code (outside
stellopt).

For input, NEO requires a Boozer coordinate file boozmn\_XXXX.nc
generated with [BOOZ\_XFORM](BOOZ_XFORM). For this tutorial we will use
the NCSX Boozer coordinate file created in the
[BOOZ\_XFORM tutorial](Boozer Transformation for NCSX-Like Configuration).
For convenience, you can also obtain the file here:
[examples/boozmn\_ncsx\_c09r00\_free.nc](examples/boozmn_ncsx_c09r00_free.nc)

------------------------------------------------------------------------

\> 1.\_\_**Edit the input text file if needed.**\_\_ \> When NEO is run
from within [STELLOPT](STELLOPT), the NEO input parameters are provided
using a fortran namelist. However, when NEO is run as a standalone code,
the input parameters are provided using a text file that is not in
fortran namelist format. This input file can be given the name
neo\_in.extension, where \'extension\' is the same extension associated
with the vmec and boozmn files. This file is read by
neo\_read\_control.f90. For this tutorial, we will use the following
input file: <file:neo_in.ncsx_c09r00_free> also shown here:

    this line is ignored
    this line is ignored
    this line is ignored
    this line is ignored
    neo_out
    10 ! no_fluxs = Number of flux surfaces to analyze.
    9 19 29 39 49 59 69 79 89 99  ! fluxs_arr = Indices of flux surfaces to analyze.
    200 ! theta_n
    200 ! phi_n
    0 ! max_m_mode
    0 ! max_n_mode
    75 ! npart
    1 ! multra
    0.01 ! acc_req
    100 ! no_bins
    75 ! nstep_per
    500 ! nstep_min
    2000 ! nstep_max
    0 ! calc_nstep_max
    1 ! eout_swi
    0 ! lab_swi
    0 ! inp_swi
    2 ! ref_swi
    1 ! write_progress
    0 ! write_output_files
    0 ! spline_test
    0 ! write_integrate
    0 ! write_diagnostic
    this line is ignored
    this line is ignored
    this line is ignored
    0 ! calc_cur
    neo_cur
    0 ! npart_cur
    0 ! alpha_cur
    0 ! write_cur_inte

\> 2. \_\_**Execute the code.**\_\_ \> To run NEO, call xneo, using the
\'extension\' as an argument. For the settings used here, there will be
nothing printed to standard output. Execution should take a few tens of
seconds.

    >~/bin/xneo ncsx_c09r00_free                                               

\> 3. \_\_**Examine the output.**\_\_ \> For this example, two output
files are created: neo\_out and neolog.ncsx\_c90r00\_free. The contents
of neo\_out should resemble the following:

            9  0.8973041424E-04  0.4370353465E+00  0.3939978795E+00  0.1601924599E+01  0.1492993106E+01
           19  0.6914966630E-04  0.8181774555E+00  0.4389194859E+00  0.1638204348E+01  0.1492993106E+01
           29  0.8145391202E-04  0.1120313880E+01  0.4759618431E+00  0.1668479141E+01  0.1492993106E+01
           39  0.1141052616E-03  0.1376483635E+01  0.5108445069E+00  0.1696584674E+01  0.1492993106E+01
           49  0.1608986166E-03  0.1601209529E+01  0.5469579776E+00  0.1723818408E+01  0.1492993106E+01
           59  0.2594685067E-03  0.1802365936E+01  0.5849734961E+00  0.1754976372E+01  0.1492993106E+01
           69  0.4658703366E-03  0.1984891866E+01  0.6216997349E+00  0.1792637760E+01  0.1492993106E+01
           79  0.7215466287E-03  0.2152625037E+01  0.6516614949E+00  0.1831918755E+01  0.1492993106E+01
           89  0.1096822199E-02  0.2308678197E+01  0.6662392549E+00  0.1872206513E+01  0.1492993106E+01
           99  0.1853340819E-02  0.2455739731E+01  0.6550476182E+00  0.1913769437E+01  0.1492993106E+01

\> The meaning of the various columns in this file is discussed on
[this page](NEO). Since we had set EOUT\_SWI=1 in the input file, the
first two columns of neo\_out contain the key information: the flux
surface index and epsilon\_effective respectively.
