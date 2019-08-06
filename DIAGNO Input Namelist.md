DIAGNO Input Namelist
=====================

------------------------------------------------------------------------

The DIAGNO input namelist is stored in the file \'diagno.control\' and
specifies the location of various files and parameters.

\|\| Variable \|\| Type \|\| Size \|\| Description \|\| \|\| nu \|\|
integer \|\| 1 \|\| Number of poloidal points to use for real-space
coordinate transformation on potential surface. \|\| \|\| nv \|\|
integer \|\| 1 \|\| Number of toroidal points to use for real-space
coordinate transformation on potential surface. \|\| \|\| input\_form
\|\| character \|\| 50 \|\| Input format specifcation (\'vmec2000\' or
\'nemec93\') \|\| \|\| diagno\_coils\_file \|\| character \|\| 256 \|\|
Machine coils file. \|\| \|\| machine\_string \|\| character \|\| 50
\|\| Machine identifier \|\| \|\| vert\_field\_file \|\| character \|\|
256 \|\| Machine verticle coils file. \|\| \|\| old\_input\_form \|\|
integer \|\| 1 \|\| Set to 1 to use old input format. \|\| \|\|
old\_cnv\_form \|\| integer \|\| 1 \|\| Set to 1 to use old cnv format.
\|\| \|\| flux\_diag\_file \|\| character \|\| 256 \|\| Flux Loop
Specification File \|\| \|\| bprobes\_file \|\| character \|\| 256 \|\|
Magnetic Field Probes Specification File \|\| \|\| mir1\_file \|\|
character \|\| 256 \|\| Mirnov Coils Specification File \|\| \|\|
seg\_rog\_file \|\| character \|\| 256 \|\| Segmented Rogowski Coil
Specification File \|\| \|\| bfield\_points\_file \|\| character \|\|
256 \|\| Magnetic Field Test Points Specification File \|\| \|\|
barrow\_param\_file \|\| character \|\| 256 \|\| Barrow Parameter
Specification File \|\| \|\| fdb\_fluxloops \|\| character \|\| 256 \|\|
Flux Loops Database File \|\| \|\| fdb\_SegRogCoils \|\| character \|\|
256 \|\| Segmented Rogowski Coil Database File \|\| \|\| lwrpl\_surf
\|\| logical \|\| 1 \|\| Flag to control output of plasma surface
quantities (execution halts after output). \|\| \|\| ltrace\_progress
\|\| logical \|\| 1 \|\| Flag to control verbose output to screen. \|\|
\|\| flux\_turns \|\| real(12,100) \|\| 100 \|\| Number of turns for
each flux loop. \|\| \|\| lcomp\_dia \|\| logical \|\| 1 \|\| Switch for
diamagnetic compensation of flux loops \|\| \|\| dia\_comp\_coefs \|\|
real(12,100) \|\| 100 \|\| Coefficients for diamagnetic compensation of
flux loops \|\| \|\| seg\_rog\_turns \|\| real(12,100) \|\| 100 \|\|
Number of turns for each segment of the segmented Rogowski coils \|\|
\|\| lflux\_comb \|\| logical \|\| 1 \|\| Switch for flux loop
combination \|\| \|\| flux\_comb\_coefs \|\| real(12,100) \|\| 100 \|\|
Coefficients for flux loop combination \|\| \|\| \<span style=\"color:
\#ff0000;\"\>lvacuum\</span\> \|\| \<span style=\"color:
\#ff0000;\"\>logical\</span\> \|\| \<span style=\"color:
\#ff0000;\"\>1\</span\> \|\| \<span style=\"color: \#ff0000;\"\>Flag to
control direct Biot-Savart calculation for vacuum data\</span\> \|\|
\|\| \<span style=\"color: \#ff0000;\"\>lcalc\_induction\</span\> \|\|
\<span style=\"color: \#ff0000;\"\>logical\</span\> \|\| \<span
style=\"color: \#ff0000;\"\>1\</span\> \|\| \<span style=\"color:
\#ff0000;\"\>Flag to control calculation of flux response file\</span\>
\|\| \|\| \<span style=\"color: \#ff0000;\"\>induction\_file\</span\>
\|\| \<span style=\"color: \#ff0000;\"\>character\</span\> \|\| \<span
style=\"color: \#ff0000;\"\>256\</span\> \|\| \<span style=\"color:
\#ff0000;\"\>Flux response file\</span\> \|\| \|\| \<span style=\"color:
\#ff0000;\"\>lnovmecb\</span\> \|\| \<span style=\"color:
\#ff0000;\"\>logical\</span\> \|\| \<span style=\"color: \#ff0000;\"\>1
\|\| \<span style=\"color: \#ff0000;\"\>Do not use VMEC B-Field
(calculate from Biot-Savart, old way)\</span\> \|\| Namelist variables
in \<span style=\"color: \#ff0000;\"\>red\</span\> are in development by
S. Lazerson.
