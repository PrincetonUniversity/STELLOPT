VMEC2XGC
========

[toc](toc) The VMEC2XGC code generates an XGC grid from a VMEC output
file.

------------------------------------------------------------------------

### Theory\[\[\#Theory\]\]

The XGC code requires a computational grid which is field aligned and
the VMEC coordinates system can easily be mapped to PEST straight field
line coordinates. The code begins by adjusting the radial resolution of
the grid so as to guarantee a minimum grid resolution (usually dictated
by the gyro-radius). A grid is then generated in the phi=0 plane which
is equidistant in the PEST poloidal angle. A reverse lookup is used to
map a given rho and u\* into the VMEC poloidal angle. Specifically the
following equation is inverted: [math](math) \\theta\^\* = \\theta +
\\lambda\\left(\\rho,\\theta,\\zeta\\right) [math](math) Here lambda is
the function which allows a transformation from PEST the VMEC poloidal
angle. In both systems the toroidal angle is the same. Once the VMEC
poloidal angle is found, the code inverts the coordinates, returning the
equilibrium quantities. The PEST poloidal angle is advanced using:
[math](math) \\theta\^\* = \\theta\^\* + \\iota\\zeta [math](math) where
iota is the rotational transform from VMEC.

------------------------------------------------------------------------

### Compilation\[\[\#Compilation\]\]

The code is compiled as part of the STELLOPT family of codes.

------------------------------------------------------------------------

### Input Data Format\[\[\#Input\]\]

This code takes a VMEC \'wout\' file as input with a few options for
controlling the grid resolution.

------------------------------------------------------------------------

### Execution\[\[\#Execution\]\]

The VMEC2XGC code is controlled through a set of command-line inputs.
[code format=\"bash\"](code format="bash") XVMEC2XGC -vmec \<VMEC FILE\>
-rad \<RADIAL RESOLUTION\> -nu \<POLOIDAL GRID\> -nv \<TOROIDAL GRID\>
-help [code](code) \|\| Argument \|\| Default \|\| Description \|\| \|\|
\<span style=\"color: \#0000ff;\"\>-vmec\</span\> \|\| \<span
style=\"color: \#0000ff;\"\>NONE\</span\> \|\| \<span style=\"color:
\#0000ff;\"\>VMEC input extension\</span\> \|\| \|\| \<span
style=\"color: \#0000ff;\"\>-rad\</span\> \|\| \<span style=\"color:
\#0000ff;\"\>0.1\</span\> \|\| \<span style=\"color: \#0000ff;\"\>Max
radial grid spacing \[m\]\</span\> \|\| \|\| \<span style=\"color:
\#0000ff;\"\>-nu\</span\> \|\| \<span style=\"color:
\#0000ff;\"\>360\</span\> \|\| \<span style=\"color: \#0000ff;\"\>Number
of poloidal gridpoints\</span\> \|\| \|\| \<span style=\"color:
\#0000ff;\"\>-nv\</span\> \|\| \<span style=\"color:
\#0000ff;\"\>360\</span\> \|\| \<span style=\"color: \#0000ff;\"\>Number
of toroidal gridpoints\</span\> \|\| \|\| \<span style=\"color:
\#0000ff;\"\>-help\</span\> \|\| \<span style=\"color:
\#0000ff;\"\>NONE\</span\> \|\| \<span style=\"color: \#0000ff;\"\>Print
help message\</span\> \|\| In it\'s simplest invokation the code
requires a VMEC wout file. [code format=\"bash\"](code format="bash")
\>/bin/xvmec2xgc -vmec ncsx\_c09r00\_free -rad 0.02 -nu 90 -nv 60
VMEC2XGC Version 1.00 \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-- EQUILIBRIUM INFO
\-\-\-\-\-\-\-\-\-\-\-\-\-\-- ASPECT RATIO: 4.469 BETA: 0.041 (total)
0.561 (poloidal) 0.044 (toroidal) TORIDAL CURRENT: -0.178653181069E+06
TORIDAL FLUX: 0.497 VOLUME: 2.962 MAJOR RADIUS: 1.442 MINOR\_RADIUS:
0.323 STORED ENERGY: 0.379813788851E-02
\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\--
\-\-\-\-\-\-\-\-\-\-\-\-\-- Adapting Radial Mesh
\-\-\-\-\-\-\-\-\-\-\-\-- NRAD DL MIN\_TARG Final NRAD: 99 DL: 0.0174
\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\--
\-\-\-\-\-\-\-\-\-\-\-\-\-- Generating R/Z Grid
\-\-\-\-\-\-\-\-\-\-\-\-\-- NRAD: 99 NU: 90 NV: 60 FILENAME:
xgc\_grid.ncsx\_c09r00\_free
\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\--
[code](code)

------------------------------------------------------------------------

### Output Data Format\[\[\#Output\]\]

The grid is output in a text file called xgc\_grid.\<ext\> where \<ext\>
is the same extension as the VMEC wout file which produced it. It is a
table with the following format: [code](code) Surf Rho Theta-VMEC PHI
Theta-PEST R \[m\] Lambda Z\[m\] Br Bphi Bz Jr Jphi Jz 1
0.0000000000000000 0.0000000000000000 0.0000000000000000
0.0000000000000000 1.6078237350321833 0.0000000000000000
0.0000000000000000 -1.40216829709591508E-008 1.4873907034966531
0.15907275299008100 0.0000000000000000 0.0000000000000000
0.0000000000000000 2 1.04123277699063634E-004 0.0000000000000000
0.0000000000000000 0.0000000000000000 1.6120197762453949
0.0000000000000000 0.0000000000000000 -1.38779083742674311E-008
1.4892487706292936 0.15809891138726201 0.0000000000000000
0.0000000000000000 0.0000000000000000 3 4.16493110796254534E-004
0.0000000000000000 0.0000000000000000 0.0000000000000000
1.6131967921932790 0.0000000000000000 0.0000000000000000
-1.37653340455094756E-008 1.4899809153977717 0.15709728550029511
0.0000000000000000 0.0000000000000000 0.0000000000000000 4
9.37109556311321867E-004 0.0000000000000000 0.0000000000000000
0.0000000000000000 1.6141917288247565 0.0000000000000000
0.0000000000000000 -1.36854708548708482E-008 1.4902683030924224
0.15622866953383502 0.0000000000000000 0.0000000000000000
0.0000000000000000 5 1.66597244318501814E-003 0.0000000000000000
0.0000000000000000 0.0000000000000000 1.6151426902998058
0.0000000000000000 0.0000000000000000 -1.36256480927259319E-008
1.4902037032805191 0.15539776420788576 0.0000000000000000
0.0000000000000000 0.0000000000000000 6 2.60308222757534014E-003
0.0000000000000000 0.0000000000000000 0.0000000000000000
1.6161069466965616 0.0000000000000000 0.0000000000000000
-1.35823167257761544E-008 1.4898102306162608 0.15456465961596752
0.0000000000000000 0.0000000000000000 0.0000000000000000 7
3.74843822524528747E-003 0.0000000000000000 0.0000000000000000
0.0000000000000000 1.6171170581109500 0.0000000000000000
0.0000000000000000 -1.35545072441891477E-008 1.4890902798382992
0.15370661575907446 0.0000000000000000 0.0000000000000000
0.0000000000000000 8 5.10204127248453654E-003 0.0000000000000000
0.0000000000000000 0.0000000000000000 1.6181950702874757
0.0000000000000000 0.0000000000000000 -1.35429774360054998E-008
1.4880379084585404 0.15280822006437170 0.0000000000000000
0.0000000000000000 0.0000000000000000 9 6.66388977274007255E-003
0.0000000000000000 0.0000000000000000 0.0000000000000000
1.6193576776912779 0.0000000000000000 0.0000000000000000
-1.35489303857185815E-008 1.4866433149846929 0.15185772227496297
0.0000000000000000 0.0000000000000000 0.0000000000000000 10
8.43398600680189681E-003 0.0000000000000000 0.0000000000000000
0.0000000000000000 1.6206185270828133 0.0000000000000000
0.0000000000000000 -1.35742724274744703E-008 1.4848948095435375
0.15084535409274014 0.0000000000000000 0.0000000000000000
0.0000000000000000 [code](code)

------------------------------------------------------------------------

### Visualization\[\[\#Visualization\]\]

The output table can be easily imported into various visualization
packages for plotting.

------------------------------------------------------------------------

### Tutorials\[\[\#Tutorials\]\]

The the above execution section for an example on how to run the code.
