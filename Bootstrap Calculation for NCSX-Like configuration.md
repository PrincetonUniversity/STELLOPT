Tutorial: BOOTSJ Calculation
============================

This tutorial will walk the user through running the [BOOTSJ](BOOTSJ)
code. For this example the National Compact Stellarator Experiment
(NCSX) configuration will be used. This machine is stellarator symmetric
with a periodicity of three. To generate the necessary files see the
tutorials [VMEC Fixed Boundary Tutorial](VMEC Fixed Boundary Run) and
[Boozer Transform Tutorial](Boozer Transformation for NCSX-Like Configuration).
You will need a \'boozmn\' file to run this simulation and the VMEC
\'input\' file.

------------------------------------------------------------------------

1.  \_\_**Edit the input namelist text file.**\_\_ \> The BOOTSJ code
    will look in the input file for the BOOTIN namelist. We will make
    use of the TEMPRES variable to choose the temperatures to be
    consistent with the [VMEC](VMEC) pressure profile. Thus the density
    will be constant.

<!-- -->

    #!fortran
    &INDATA
      MGRID_FILE = 'none'
      LFREEB = F
      DELT =   9.00E-01
      TCON0 =   2.00E+00
      .
      .
      .
      RBC(1,10) = -1.08453911913102E-06     ZBS(1,10) =  8.60890353665843E-05
      RBC(2,10) =  1.04051545504927E-04     ZBS(2,10) = -2.17661420286656E-04
      RBC(3,10) = -5.21965328013036E-04     ZBS(3,10) = -2.67111216700977E-04
      RBC(4,10) = -4.95991087393098E-04     ZBS(4,10) =  2.43875640076056E-05
      RBC(5,10) = -1.94520415280627E-04     ZBS(5,10) =  1.55759001593971E-04
      RBC(6,10) = -6.94143617569942E-05     ZBS(6,10) =  4.40565098025554E-05
    /
    &BOOTIN
      MBUSE = 72  ! From Boozer Transformation
      NBUSE = 24  ! From Boozer Transformation
      ZEFF1 = 1.0
      DAMP_BS = 0.001 
      TEMPRES = 1.0 ! Coefficient Te=P^TEMPRES
      TETI    = 1.0 ! Te/Ti
      DENS0   = 1.0 ! 10^20 [m^-3]
    /

1.  **\_\_Edit the input text file\_\_** You will also need to create a
    test file with the extension of the Boozer transformation run and
    the surfaces on which you wish to compute the bootstrap current.

<!-- -->

    ncsx_c09r00_fixed
    2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 
    32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 
    60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 
    88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99

1.  \_\_**Execute the code.**\_\_ You simply need to pass input text
    filename on the command line to BOOTSJ. You will need to make sure
    you have the \'boozmn\' file in the directory as well.

<!-- -->

    ~/bin/xbootsj bootsj_in.ncsx_c09r00_fixed
     Start BOOTSJ: Version 7.01 
     mbuse =           72 nbuse =           24 nthetah =          192  nzetah =           96
     Total bootstrap current =   0.0162554 MA
     Finished BOOTSJ, time =   119.740  sec

1.  \_\_**Examine the output.**\_\_ This simulation outputs three files:
    answers.\<ext\>, answers\_plot.\<ext\> and jBbs.\<ext\>. These files
    contain various tables in text form relating to the bootstrap
    current and run parameters.