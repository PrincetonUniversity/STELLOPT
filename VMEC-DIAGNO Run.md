Tutorial: VMEC-DIAGNO Run
=========================

This tutorial will guide the user through setup and execution of the
DIAGNO routine given a VMEC output. DIAGNO is more commonly used as a
module in the STELLOPT optimization routine and it is rare that the user
would actually want to run DIAGNO as a stand-alone application. This
tutorial focuses more on the proper setup of the input files so the user
has an idea of what is required when running STELLOPT in reconstruction
mode.

------------------------------------------------------------------------

1.  **\_\_Execute VMEC with the DIAGNO flag.\_\_** \> The VMEC code will
    produce a DIAGNO input file if the LDIAGNO namelist parameter is set
    to T (true) and VMEC is run in free boundary mode. The first few
    lines of an example VMEC input file is shown: \>
    [code format=\"fortran\"](code format="fortran") &INDATA MGRID\_FILE
    = \'mgrid.sample\' LFREEB = T LDIAGNO = T DELT = 9.00E-01 TCON0 =
    2.00E+00 NFP = 3 NCURR = 1 MPOL = 11 NTOR = 6 NZETA = 32 . . .
    [code](code) We have also supplied an mgrid file for the free
    boundary run. VMEC will then create a DIAGNO input file
    diagno\_in.ext (where ext is the suffix of your VMEC input
    namelist). An example of this file can be found
    here(<file:diagno_in.val080>). This is a text file containing a
    specification of the outer most flux surface, the potential on that
    surface, the magnetic axis, the current on that axis, and the
    current in each external coil group for that run.
2.  **\_\_Create the magnetic diagnostics files.\_\_** \> The DIAGNO
    routine requires that the user specify the magnetic diagnostics of
    interest for a given run. It need not require a specification for
    all diagnostics, only those for which the user wishes output. The
    specifications fall into two categories: point measurements and
    integrated loops. For point measurements we simply need to specify
    the location of a measurement relative to the center of the device
    (x,y,z) and the orientation of the measurement (theta,zeta). The
    bpoints file only requires the position be specified as it only
    produces log-file output of the field at a point (used for testing).
    Integrated Loops require the user specify the loop as a series of
    points along the loop (x,y,z). This is similar to the specification
    of magnetic coils only the current carried is not used.
3.  \_\_**Edit the diagno.control file.**\_\_ \> The input namelist for
    the DIAGNO routine is stored in the \'diagno.control\' file. The
    [DIAGNO input namelist](DIAGNO input namelist) controls how DIAGNO
    runs. If a diagnostic is not specified in this file then it is not
    evaluated by DIAGNO. Here is an example file: \>
    [code format=\"fortran\"](code format="fortran") &DIAGNO\_IN nu =
    360 nv = 36 input\_form = \'vmec2000\' machine\_string = \'EXMAPLE\'
    diagno\_coils\_file = \'coils.example\' flux\_diag\_file =
    \'test\_diagno.flux\_file\' flux\_turns = 1 -1 1 -1 1 -1 1 -1 -1 1 1
    -1 1 1 1 1 1 1 1 -1 1 -1 1 -1 lcomp\_dia = .false. lflux\_comb =
    .false. ltrace\_progress = .true. lwrpl\_surf = .false. / &END
    [code](code)
4.  \_\_**Execute DIAGNO.**\_\_ \> To simplify execution of the code,
    the VMEC compilation scripts create a directory called \'bin\' in
    your home () directory. Symbolic links are then placed there
    pointing to each of the compiled codes in their respective
    \'Vrelease\' subdirectories. In practice, the screen output from
    DIAGNO should be redirected to a log file and put in the background
    (\>& log\_diagno.example &) but for this tutorial we will simply
    instruct the user to run the code so they can see how DIAGNO
    executes. This is done by passing the suffix of the VMEC output file
    to the DIAGNO code through the command line. The diagno.control file
    must be in the same directory. \> [code](code) \> /bin/xdiagno
    example New input form! vmec2000 Machine to be handled: EXAMPLE
    Number of field periods: 5 Coil currents LHD IHCI/A IHCM/A IHCO/A
    IPIV/A IPIS/A IPOV/A 4.8325E+03 3.1647E+03 1.9311E+03 -5.3180E+02
    1.0000E+00 5.8200E+02 tor.curr/A 3.0655E-18 nu nv nuv nvp nuvp 360
    36 12960 360 129600 Coils of LHD requested! First Pass through coils
    Number of elements: 752324 Allocating arrays Number of Currents: 6
    Reading Coils nels: 108300 k: 108301 curre: 4832.500000000000 nels:
    108300 k: 216601 curre: 3164.700000000000 nels: 108300 k: 324901
    curre: 1931.100000000000 nels: 103968 k: 428869 curre:
    -531.8000000000000 nels: 150176 k: 579045 curre: 1.000000000000000
    nels: 173280 k: 752325 curre: 582.0000000000000 Calling
    initialize\_surface\_values! Calculating Trig. Functions Calculating
    potential in real space Calculating Magnetic Axis rmin at phi= 0
    degrees= 4.0311E+00 rmax at phi= 0 degrees= 4.1851E+00 zmax at phi=
    0 degrees= 4.9736E-01 rmin at phi=36 degrees= 4.1475E+00 rmax at
    phi=36 degrees= 4.4908E+00 zmax at phi=36 degrees= 3.4088E-01 rav=
    2.22625E-01m Calculating components of outward normal Magnetic Field
    due to Coils phi=0: x y z Bx By Bz 4.17850E+00 0.00000E+00
    0.00000E+00 -2.27380E-18 7.88501E-01 3.38170E-01 4.17741E+00
    0.00000E+00 2.98786E-01 1.59822E-01 7.24842E-01 3.33743E-01
    4.18369E+00 0.00000E+00 4.84250E-01 2.48969E-01 6.25666E-01
    3.23898E-01 4.10653E+00 0.00000E+00 4.65334E-01 2.14802E-01
    6.38228E-01 2.93470E-01 4.03212E+00 0.00000E+00 2.99086E-01
    1.29103E-01 7.13479E-01 2.66851E-01 4.04066E+00 0.00000E+00
    6.08676E-17 2.01833E-17 7.65942E-01 2.70954E-01 4.03212E+00
    0.00000E+00 -2.99086E-01 -1.29103E-01 7.13479E-01 2.66851E-01
    4.10653E+00 0.00000E+00 -4.65334E-01 -2.14802E-01 6.38228E-01
    2.93470E-01 4.18369E+00 0.00000E+00 -4.84250E-01 -2.48969E-01
    6.25666E-01 3.23898E-01 4.17741E+00 0.00000E+00 -2.98786E-01
    -1.59822E-01 7.24842E-01 3.33743E-01 Calculating Vec. Potential due
    to plasma. \--Adding Coil Contribution \--Calling becoil quite
    often! Finishing calculation and extending over whole plasma
    boundary! Freeing local allocations in surface values! Initializing
    surface values finished! Before diagnostics! 24 flux loops FluxLoop(
    1 ) Segments: 16 Toroidal Flag: 0 Diamag. Flag: 0 Name:VSL0011
    FluxLoop( 2 ) Segments: 16 Toroidal Flag: 0 Diamag. Flag: 0
    Name:VSL0012 FluxLoop( 3 ) Segments: 16 Toroidal Flag: 0 Diamag.
    Flag: 0 Name:VSL0013 FluxLoop( 4 ) Segments: 16 Toroidal Flag: 0
    Diamag. Flag: 0 Name:VSL0021 FluxLoop( 5 ) Segments: 16 Toroidal
    Flag: 0 Diamag. Flag: 0 Name:VSL0022 FluxLoop( 6 ) Segments: 16
    Toroidal Flag: 0 Diamag. Flag: 0 Name:VSL0023 FluxLoop( 7 )
    Segments: 16 Toroidal Flag: 0 Diamag. Flag: 0 Name:VSL0031 FluxLoop(
    8 ) Segments: 16 Toroidal Flag: 0 Diamag. Flag: 0 Name:VSL0032
    FluxLoop( 9 ) Segments: 16 Toroidal Flag: 0 Diamag. Flag: 0
    Name:VSL0033 FluxLoop( 10 ) Segments: 16 Toroidal Flag: 0 Diamag.
    Flag: 0 Name:VSL0041 FluxLoop( 11 ) Segments: 16 Toroidal Flag: 0
    Diamag. Flag: 0 Name:VSL0042 FluxLoop( 12 ) Segments: 16 Toroidal
    Flag: 0 Diamag. Flag: 0 Name:VSL0043 FluxLoop( 13 ) Segments: 16
    Toroidal Flag: 0 Diamag. Flag: 0 Name:VSL0051 FluxLoop( 14 )
    Segments: 16 Toroidal Flag: 0 Diamag. Flag: 0 Name:VSL0052 FluxLoop(
    15 ) Segments: 16 Toroidal Flag: 0 Diamag. Flag: 0 Name:VSL0053
    FluxLoop( 16 ) Segments: 16 Toroidal Flag: 0 Diamag. Flag: 0
    Name:VSL0061 FluxLoop( 17 ) Segments: 16 Toroidal Flag: 0 Diamag.
    Flag: 0 Name:VSL0062 FluxLoop( 18 ) Segments: 16 Toroidal Flag: 0
    Diamag. Flag: 0 Name:VSL0063 FluxLoop( 19 ) Segments: 16 Toroidal
    Flag: 0 Diamag. Flag: 0 Name:VSL0071 FluxLoop( 20 ) Segments: 16
    Toroidal Flag: 0 Diamag. Flag: 0 Name:VSL0072 FluxLoop( 21 )
    Segments: 16 Toroidal Flag: 0 Diamag. Flag: 0 Name:VSL0073 FluxLoop(
    22 ) Segments: 16 Toroidal Flag: 0 Diamag. Flag: 0 Name:VSL0081
    FluxLoop( 23 ) Segments: 16 Toroidal Flag: 0 Diamag. Flag: 0
    Name:VSL0082 FluxLoop( 24 ) Segments: 16 Toroidal Flag: 0 Diamag.
    Flag: 0 Name:VSL0083 Maximum number of segments+1 = 17 Reading Flux
    Loops Device:LHD converting mm-\>m Performing line integral using
    closed simpson formula coil number 1 VSL0011 nseg= 16 flux( 1)=
    1.7625074025E-05 Wb Flux normalized to phiedge/2 : 3.52501E-04

coil number 2 VSL0012 nseg= 16 flux( 2)= 1.1103216535E-04 Wb Flux
normalized to phiedge/2 : 2.22064E-03

coil number 3 VSL0013 nseg= 16 flux( 3)= 5.4062372699E-04 Wb Flux
normalized to phiedge/2 : 1.08125E-02

coil number 4 VSL0021 nseg= 16 flux( 4)= 2.1563993016E-03 Wb Flux
normalized to phiedge/2 : 4.31280E-02

coil number 5 VSL0022 nseg= 16 flux( 5)= 4.6510679098E-03 Wb Flux
normalized to phiedge/2 : 9.30214E-02

coil number 6 VSL0023 nseg= 16 flux( 6)= 2.7574302570E-03 Wb Flux
normalized to phiedge/2 : 5.51486E-02

coil number 7 VSL0031 nseg= 16 flux( 7)= -2.7573435220E-03 Wb Flux
normalized to phiedge/2 : -5.51469E-02

coil number 8 VSL0032 nseg= 16 flux( 8)= -4.6810560303E-03 Wb Flux
normalized to phiedge/2 : -9.36211E-02

coil number 9 VSL0033 nseg= 16 flux( 9)= -2.1562598803E-03 Wb Flux
normalized to phiedge/2 : -4.31252E-02

coil number 10 VSL0041 nseg= 16 flux(10)= 2.7570436880E-03 Wb Flux
normalized to phiedge/2 : 5.51409E-02

coil number 11 VSL0042 nseg= 16 flux(11)= 4.6815182993E-03 Wb Flux
normalized to phiedge/2 : 9.36304E-02

coil number 12 VSL0043 nseg= 16 flux(12)= 2.1560640851E-03 Wb Flux
normalized to phiedge/2 : 4.31213E-02

coil number 13 VSL0051 nseg= 16 flux(13)= -5.4064632306E-04 Wb Flux
normalized to phiedge/2 : -1.08129E-02

coil number 14 VSL0052 nseg= 16 flux(14)= -1.1101311193E-04 Wb Flux
normalized to phiedge/2 : -2.22026E-03

coil number 15 VSL0053 nseg= 16 flux(15)= -1.7617257528E-05 Wb Flux
normalized to phiedge/2 : -3.52345E-04

coil number 16 VSL0061 nseg= 16 flux(16)= 5.4062372699E-04 Wb Flux
normalized to phiedge/2 : 1.08125E-02

coil number 17 VSL0062 nseg= 16 flux(17)= 1.1103216535E-04 Wb Flux
normalized to phiedge/2 : 2.22064E-03

coil number 18 VSL0063 nseg= 16 flux(18)= 1.7625074025E-05 Wb Flux
normalized to phiedge/2 : 3.52501E-04

coil number 19 VSL0071 nseg= 16 flux(19)= -2.7573435220E-03 Wb Flux
normalized to phiedge/2 : -5.51469E-02

coil number 20 VSL0072 nseg= 16 flux(20)= -4.6810560303E-03 Wb Flux
normalized to phiedge/2 : -9.36211E-02

coil number 21 VSL0073 nseg= 16 flux(21)= -2.1562598803E-03 Wb Flux
normalized to phiedge/2 : -4.31252E-02

coil number 22 VSL0081 nseg= 16 flux(22)= -5.4060001228E-04 Wb Flux
normalized to phiedge/2 : -1.08120E-02

coil number 23 VSL0082 nseg= 16 flux(23)= -1.1100662611E-04 Wb Flux
normalized to phiedge/2 : -2.22013E-03

coil number 24 VSL0083 nseg= 16 flux(24)= -1.7948603750E-05 Wb Flux
normalized to phiedge/2 : -3.58972E-04 [code](code)

1.  \_\_**Examine the output.**\_\_ \> For each diagnostic type (the
    execption being the bfield\_points) the code will produce a
    diagnostic output file. These value should be compared against the
    measured values. Here the STELLOPT routine iteratively modifies the
    VMEC input parameters until a good fit between simulated and
    measured values is found.
