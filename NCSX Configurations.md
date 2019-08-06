NCSX Configurations
===================

Equilibria were generated using the coils file c07r00 posted in
<ftp://ftp.pppl.gov/pub/abrooks/NCSX/NCSX_COIL_DATA/li383_1.4m/coilsets>
. Plasma-to-FW (first wall) distances are calculated using fw030324.vmec
posted in
<ftp://ftp.pppl.gov/pub/abrooks/NCSX/NCSX_COIL_DATA/li383_1.4m/FirstWall>
Note: Coil Currents quoted below are STELLOPT EXTCUR values. EXTCUR(1-3)
correspond to the independent modular coils; EXTCUR(4-5) correspond to
the OH coils; EXTCUR(6-9) correspond to the PF coils EXTCUR(10)
corresponds to the TF coil To calculate Amp-turns for each coil group,
multiply by the number of turns in the appropriate group. The turns
corresponding to the 10 coil groups are NTURNS(1-3) = 1; NTURNS(4-5) =
4, 4; NTURNS(6-9) = 4, 6, 2, 1; NTURNS(10) = 1.

------------------------------------------------------------------------

Overview
--------

\|\| Configuration \|\| Modular 1 \|\| Modular 2 \|\| Modular3 \|\| PF1
\|\| PF2 \|\| PF3 \|\| PF4 \|\| PF5 \|\| PF6 \|\| TF \|\| CURTOR \|\|
BETA \|\| \|\| S1 \|\| 7.5538E+05 \|\| 7.1836E+05 \|\| 6.3899E+05 \|\| 0
\|\| 0 \|\| -1.9920E+04 \|\| 4.0100E+05 \|\| -5.3654E+04 \|\| 9.2250E+03
\|\| -4.3070E+04 \|\| 0.0 \|\| 0.0 \|\| \|\| S2a \|\| 7.0081E+05 \|\|
6.4982E+05 \|\| 5.8390E+05 \|\| 0 \|\| 0 \|\| -6.1276E+04 \|\|
-7.4436E+05 \|\| 2.6111E+05 \|\| -1.2157E+05 \|\| 1.6501E+04 \|\|
-1.786E+05 \|\| 1.0 \|\| \|\| S2b \|\| 7.0010E+05 \|\| 6.4954E+05 \|\|
5.7806E+05 \|\| 0 \|\| 0 \|\| -6.0047E+04 \|\| -7.7195E+05 \|\|
2.1485E+05 \|\| -1.1309E+05 \|\| 1.9448E+04 \|\| -1.600E+05 \|\| 0.0
\|\| \|\| S3 \|\| 6.5227E+05 \|\| 6.5187E+05 \|\| 5.3774E+05 \|\| 0 \|\|
0 \|\| 2.8095E+04 \|\| -5.4805E+04 \|\| 3.0123E+04 \|\| 9.4241E+04 \|\|
4.5514E+04 \|\| -1.786E+05 \|\| 4.1 \|\| \|\| S3-mod \|\| 6.5858E+05
\|\| 6.5440E+05 \|\| 5.4252E+05 \|\| 0 \|\| 0 \|\| 2.6137E+04 \|\|
-5.8987E+04 \|\| 2.7881E+04 \|\| 8.9996E+04 \|\| 4.5309E+04 \|\|
-1.786E+05 \|\| 4.1 \|\| \|\| S4 \|\| 6.8234E+05 \|\| 6.3959E+05 \|\|
5.7791E+05 \|\| 0 \|\| 0 \|\| -8.3748E+04 \|\| -6.0941E+05 \|\|
3.1129E+05 \|\| -1.8187E+05 \|\| 2.6233E+04 \|\| -3.200E+05 \|\| 0.0
\|\| \|\| S5 \|\| 7.0010E+05 \|\| 6.4954E+05 \|\| 5.7806E+05 \|\| 0 \|\|
0 \|\| -9.9202E+04 \|\| -3.4874E+05 \|\| 1.6463E+05 \|\| -1.9729E+05
\|\| 1.9448E+04 \|\| 0.0 \|\| 0.0 \|\|

------------------------------------------------------------------------

S1: Vac. Configuration with iota \> 0.5 everywhere in plasma
------------------------------------------------------------

\|\| TOROIDAL CURRENT\* \|\| 4.47E-16 \[MA\] \|\| \|\| R-BTOR(s=1) \|\|
2.38E+00 \|\| \|\| Average Beta \|\| 0.000E+00 \|\| \|\| Aspect Ratio
\|\| 4.682891 \|\| \|\| Plasma Volume \|\| 2.718500 \[M**3\] \|\| \|\|
Major Radius \|\| 1.445471 \[M\] (from Volume and Cross Section) \|\|
\|\| Minor Radius \|\| 0.308671 \[M\] (from Cross Section) \|\| \|\|
Waist (v = 0) \|\| 0.238598 \[M\] in R \|\| \|\| Full Height(v = 0) \|\|
1.253904 \[M\] \|\| \|\| Waist (v = pi) in R \|\| 0.740767 \[M\] \|\|
\|\| Full Height(v = pi) \|\| 0.510073 \[M\] \|\| \|\|
iota(s=.00,.50,1.00) \|\| 5.297E-01 5.374E-01 5.677E-01 \|\| \|\|
eps32(js=16,26,40) \|\| 0.10146391E-03, 0.36778945E-03, 0.16842455E-02
\|\| \|\| Min Plasma-VV Sep \|\| 2.566E-02 \|\| This was generated from
a C06R00 equilibrium generated for vacuum divertor studies of
Georgevskiy and Rudakov in Jan 03. The modulars coils for C06R00 are the
same as C07R00 so the same currents could be used. The PF coils had
changed, so the PF fields produced by the C06R00 configuration were fit
to the fields produced by appropriate C07R00 coil currents (by A. Brooks
using fitpf2.f). No further optimization was done. STELLOPT output files
stored in**

<ftp://ftp.pppl.gov/pub/pomphrey/M50.c07r00/General.Repository/m50.c07r00_s1gen.01.zip>

------------------------------------------------------------------------

S2 (full current, zero beta)
----------------------------

Using RBt = 2.38 T-m and the full plasma current (-178.6 kA), I was
unsuccessful in producing a kink stable plasma at b = 0 %. The strategy
I adopted was to produce a sequence of kink and ballooning stable states
with good FW-plasma separation as b was decreased from the S3 value of
4%. However, a strong n=1 kink instability was found to develop for b in
the range \< 0.75% at the full current preventing access to a true S2
state. The region of instability seems to be fairly localized in Ip-b
space, since a kink stable configuration was found at a slightly lower
plasma current (-160.0 kA) at b = 0%. Below, I will list two surrogate
S2 equilibria which the engineering analysis can use to develop coil
current waveforms. The first, S2a, is the low b full current equilibrium
(Ip = -178.6 kA, b = 1.0%), and the second, S2b, is the b = 0%
equilibrium with 10% reduced current (Ip = -160.0 kA, b =0%).

S2a (full current, zero beta)
-----------------------------

\|\| \*TOROIDAL CURRENT\* \|\| -1.79E-01 \[MA\] \|\| \|\|
\*R-BTOR(s=1)\* \|\| 2.38E+00 \|\| \|\| \*Average Beta\* \|\| 9.719E-03
\|\| \|\| \*Aspect Ratio\* \|\| 4.486785 \|\| \|\| \*Plasma Volume\*
\|\| 2.887553 \[M**3\] \|\| \|\| \*Major Radius\* \|\| 1.433366 \[M\]
(from Volume and Cross Section) \|\| \|\| \*Minor Radius\* \|\| 0.319464
\[M\] (from Cross Section) \|\| \|\| \*Waist (v = 0)\* \|\| 0.241542
\[M\] in R \|\| \|\| \*Full Height(v = 0)\* \|\| 1.340239 \[M\] \|\|
\|\| \*Waist (v = pi) in R\* \|\| 0.651921 \[M\] \|\| \|\| \*Full
Height(v = pi)\* \|\| 0.578833 \[M\] \|\| \|\| \*iota(s=.00,.50,1.00)\*
\|\| 3.499E-01 5.263E-01 6.587E-01 \|\| \|\| \*eps32(js=16,26,40)\* \|\|
0.44228472E-04, 0.17115496E-03, 0.77940978E-03 \|\| \|\| \*Min Plasma-VV
Sep\* \|\| 2.669E-02 \|\| Comments: Kink and ballooning stable. STELLOPT
output files stored in
<ftp://ftp.pppl.gov/pub/pomphrey/M50.c07r00/General.Repository/m50.c07r00_s2gen.05.zip>
==S2b: Ip = -160.0 kA, b = 0== \|\| \*TOROIDAL CURRENT\* \|\| -1.60E-01
\[MA\] \|\| \|\| \*R-BTOR(s=1)\* \|\| 2.38E+00 \|\| \|\| \*Average
Beta\* \|\| 0.000E+00 \|\| \|\| \*Aspect Ratio\* \|\| 4.642714 \|\| \|\|
\*Plasma Volume\* \|\| 2.805697 \[M**3\] \|\| \|\| \*Major Radius\* \|\|
1.452396 \[M\] (from Volume and Cross Section) \|\| \|\| \*Minor
Radius\* \|\| 0.312833 \[M\] (from Cross Section) \|\| \|\| \*Waist (v =
0)\* \|\| 0.240214 \[M\] in R \|\| \|\| \*Full Height(v = 0)\* \|\|
1.323896 \[M\] \|\| \|\| \*Waist (v = pi) in R\* \|\| 0.646621 \[M\]
\|\| \|\| \*Full Height(v = pi)\* \|\| 0.574251 \[M\] \|\| \|\|
\*iota(s=.00,.50,1.00)\* \|\| 3.342E-01 5.155E-01 6.568E-01 \|\| \|\|
\*eps32(js=16,26,40)\* \|\| 0.66174144E-04, 0.24348984E-03,
0.10503225E-02 \|\| \|\| \*Min Plasma-VV Sep\* \|\| 1.853E-02 \|\|
Comments: Kink stable. STELLOPT output files stored in
<ftp://ftp.pppl.gov/pub/pomphrey/M50.c07r00/General.Repository/m50.c07r00_s2gen.07.zip>

------------------------------------------------------------------------

S3: - full current (-178.6 kA), full beta (4.1%) RBt = 2.38 T-m (1.7 T at 1.4 m)
--------------------------------------------------------------------------------

\|\| \*TOROIDAL CURRENT\* \|\| -1.79E-01 \[MA\] \|\| \|\|
\*R-BTOR(s=1)\* \|\| 2.37E+00 \|\| \|\| \*Average Beta\* \|\| 4.081E-02
\|\| \|\| \*Aspect Ratio\* \|\| 4.466509 \|\| \|\| \*Plasma Volume\*
\|\| 2.959223 \[M**3\] \|\| \|\| \*Major Radius\* \|\| 1.440771 \[M\]
(from Volume and Cross Section) \|\| \|\| \*Minor Radius\* \|\| 0.322572
\[M\] (from Cross Section) \|\| \|\| \*Waist (v = 0)\* \|\| 0.268353
\[M\] in R \|\| \|\| \*Full Height(v = 0)\* \|\| 1.257869 \[M\] \|\|
\|\| \*Waist (v = pi) in R\* \|\| 0.719456 \[M\] \|\| \|\| \*Full
Height(v = pi)\* \|\| 0.526709 \[M\] \|\| \|\| \*iota(s=.00,.50,1.00)\*
\|\| 3.506E-01 5.522E-01 6.548E-01 \|\| \|\| \*eps32(js=16,26,40)\* \|\|
0.87389060E-04, 0.18022688E-03, 0.77031260E-03 \|\| \|\| \*Min Plasma-VV
Sep\* \|\| 1.889E-02 \|\| Comments: Kink stable to n=1 and n=0 families;
Ballooning unstable on surfaces 39-44 for theta = 0, zeta = 0
(apparently stable to finite high n ballooning according to Mike Z,
relating calculations of G. Fu). STELLOPT output files for S3 in
<ftp://ftp.pppl.gov/pub/pomphrey/M50.c07r00/Ref.Scenarios/m50.c07r00_s3.zip>
==S3-modified: (STELLOPT used to increase FW-plasma separation)== \|\|
\*TOROIDAL CURRENT\* \|\| -1.79E-01 \[MA\] \|\| \|\| \*R-BTOR(s=1)\*
\|\| 2.39E+00 \|\| \|\| \*Average Beta\* \|\| 4.047E-02 \|\| \|\|
\*Aspect Ratio\* \|\| 4.530665 \|\| \|\| \*Plasma Volume\* \|\| 2.909880
\[M**3\] \|\| \|\| \*Major Radius\* \|\| 1.446405 \[M\] (from Volume and
Cross Section) \|\| \|\| \*Minor Radius\* \|\| 0.319248 \[M\] (from
Cross Section) \|\| \|\| \*Waist (v = 0)\* \|\| 0.267560 \[M\] in R \|\|
\|\| \*Full Height(v = 0)\* \|\| 1.245501 \[M\] \|\| \|\| \*Waist (v =
pi) in R\* \|\| 0.717620 \[M\] \|\| \|\| \*Full Height(v = pi)\* \|\|
0.522160 \[M\] \|\| \|\| \*iota(s=.00,.50,1.00)\* \|\| 3.510E-01
5.588E-01 6.668E-01 \|\| \|\| \*eps32(js=16,26,40)\* \|\|
0.91137024E-04, 0.19314765E-03, 0.80924301E-03 \|\| \|\| \*Min Plasma-VV
Sep\* \|\| 2.669E-02 \|\| Comments: Stability properties same as
reference S3 above and eps32 comparable. However, compared with
reference S3 above, the Min Plasma-FW distance has increased to from
1.89 cm to 2.70 cm. STELLOPT output files for S3-modified in
<ftp://ftp.pppl.gov/pub/pomphrey/M50.c07r00/Ref.Scenarios/m50.c07r00_s3.02.zip>

------------------------------------------------------------------------

S4: High Current Ohmic (Ip = -320kA BT = 1.7T at R = 1.4 m; iota \~ 1.0)
------------------------------------------------------------------------

\|\| \*TOROIDAL CURRENT\* \|\| -3.20E-01 \[MA\] \|\| \|\|
\*R-BTOR(s=1)\* \|\| 2.38E+00 \|\| \|\| \*Average Beta\* \|\| 0.000E+00
\|\| \|\| \*Aspect Ratio\* \|\| 4.569327 \|\| \|\| \*Plasma Volume\*
\|\| 2.853251 \[M**3\] \|\| \|\| \*Major Radius\* \|\| 1.445123 \[M\]
(from Volume and Cross Section) \|\| \|\| \*Minor Radius\* \|\| 0.316266
\[M\] (from Cross Section) \|\| \|\| \*Waist (v = 0)\* \|\| 0.230570
\[M\] in R \|\| \|\| \*Full Height(v = 0)\* \|\| 1.231942 \[M\] \|\|
\|\| \*Waist (v = pi) in R\* \|\| 0.815659 \[M\] \|\| \|\| \*Full
Height(v = pi)\* \|\| 0.542858 \[M\] \|\| \|\| \*iota(s=.00,.50,1.00)\*
\|\| 3.685E-01 6.634E-01 9.607E-01 \|\| \|\| \*eps32(js=16,26,40)\* \|\|
0.63558447E-04, 0.29107842E-03, 0.12429393E-02 \|\| \|\| \*Min Plasma-VV
Sep\* \|\| 5.592E-03 \|\|**

Comments: The plasma current here is 320kA, not the 350kA we had at the
CDR. The toroidal field is 1.7T, not 1.8T. However 1.8/1.7 x 320 = 350
so the scaled field are equivalent. The edge iota is 0.96 (see above
Table), and the max iota is found to be 0.99 (at s = 0.90). For this S4
plasma shape, we find the vacuum contribution to iota at the axis and
edge are is4vac(0) = 0.5230 and is4vac(1) = 0.6329 ==\> the plasma
contributions to iota are is4plas(0) = --0.1545 and is4plas(1) =
+0.3278. These values can be compared with the corresponding
contributions for the s3 reference plasma, namely is3vac(0) = 0.3910,
is3vac(1) = 0.4857, is3plas(0) = -0.0404, is3plas(1) = +0.1691. So the
s4 state has essentially doubled the plasma contribution to iota
relative to the S3 state. Output files for the S4 state can be found in
<ftp://ftp.pppl.gov/pub/pomphrey/M50.c07r00/General.Repository/m50.c07r00_s4gen.104.zip>

------------------------------------------------------------------------

S5: Vac. configuration with same toroidal (modular + TF) coil currents as S2
----------------------------------------------------------------------------

\|\| \*TOROIDAL CURRENT\* \|\| 1.80E-16 \[MA\] \|\| \|\| \*R-BTOR(s=1)\*
\|\| 2.38E+00 \|\| \|\| \*Average Beta\* \|\| 0.000E+00 \|\| \|\|
\*Aspect Ratio\* \|\| 4.670057 \|\| \|\| \*Plasma Volume\* \|\| 2.841449
\[M**3\] \|\| \|\| \*Major Radius\* \|\| 1.464260 \[M\] (from Volume and
Cross Section) \|\| \|\| \*Minor Radius\* \|\| 0.313542 \[M\] (from
Cross Section) \|\| \|\| \*Waist (v = 0)\* \|\| 0.246149 \[M\] in R \|\|
\|\| \*Full Height(v = 0)\* \|\| 1.246927 \[M\] \|\| \|\| \*Waist (v =
pi) in R\* \|\| 0.666425 \[M\] \|\| \|\| \*Full Height(v = pi)\* \|\|
0.557373 \[M\] \|\| \|\| \*iota(s=.00,.50,1.00)\* \|\| 3.795E-01
3.914E-01 4.201E-01 \|\| \|\| \*eps32(js=16,26,40)\* \|\|
0.20579201E-03, 0.68248127E-03, 0.22441680E-02 \|\| \|\| \*Min Plasma-VV
Sep\* \|\| 2.693E-02 \|\|**

Comments: STELLOPT output files stored in
<ftp://ftp.pppl.gov/pub/pomphrey/M50.c07r00/General.Repository/m50.c07r00_s5gen.03.zip>

------------------------------------------------------------------------

Flexibility Studies
-------------------

The following files reference the Flexibility Study done for NCSX ().
<file:NCSX_Magnetic_Configuration_Flexibility_Pomphrey.pdf>

### Reference Configuration S3 (Table I & II)

<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/REF.SCENARIOS/z16.s3.1.zip>

### Coil Perturbation Study (Table IV)

\|\| Perturbed Quantitity (5%) \|\| File Location \|\| \|\| Modular 1
\|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/ROBUSTNESS/s3m1+5%25.zip>
\|\| \|\| Modular 2 \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/ROBUSTNESS/s3m2+5%25.zip>
\|\| \|\| Modular 3 \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/ROBUSTNESS/s3m3+5%25.zip>
\|\| \|\| Poloidal Field 3 \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/ROBUSTNESS/s3pf3+5%25.zip>
\|\| \|\| Poloidal Field 4 \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/ROBUSTNESS/s3pf4+5%25.zip>
\|\| \|\| Poloidal Field 5 \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/ROBUSTNESS/s3pf5+5%25.zip>
\|\| \|\| Poloidal Field 6 \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/ROBUSTNESS/s3pf6+5%25.zip>
\|\| \|\| Toroidal Field \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/ROBUSTNESS/s3tf+5%25.zip>
\|\| \|\| Plasma Current \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/ROBUSTNESS/s3ip+5%25.zip>
\|\|

### Plasma Current / Beta Scan (Table V)

\|\| \|\| \|\| \|\| Beta \|\| \|\| \|\| \|\| \|\| \|\| Ip \|\| 0.0 \|\|
1.0 \|\| 2.0 \|\| 3.0 \|\| 4.0 \|\| 5.0 \|\| 6.0 \|\| \|\| 0 \|\|
[z16.a0g0b0i0.2.zip](ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/Ip-BETA/z16.a0g0b0i0.2.zip)
\|\|
[z16.a0g0b1i0.0.zip](ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/Ip-BETA/z16.a0g0b1i0.0.zip)
\|\|
[z16.a0g0b2i0.2.zip](ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/Ip-BETA/z16.a0g0b2i0.2.zip)
\|\| Unstable \|\| Unstable \|\| XX \|\| XX \|\| \|\| 44 \|\|
[z16.a0g0b0i44.0.zip](ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/Ip-BETA/z16.a0g0b0i44.0.zip)
\|\|
[z16.a0g0b1i44.0.zip](ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/Ip-BETA/z16.a0g0b1i44.0.zip)
\|\|
[z16.a0g0b2i44.0.zip](ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/Ip-BETA/z16.a0g0b2i44.0.zip)
\|\|
[z16.a0g0b3i44.0.zip](ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/Ip-BETA/z16.a0g0b3i44.0.zip)
\|\|
[z16.a0g0b4i44.0.zip](ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/Ip-BETA/z16.a0g0b4i44.0.zip)
\|\| XX \|\| XX \|\| \|\| 87.5 \|\|
[z16.a0g0b0i875.2.zip](ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/Ip-BETA/z16.a0g0b0i875.2.zip)
\|\|
[z16.a0g0b1i875.1.zip](ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/Ip-BETA/z16.a0g0b1i875.1.zip)
\|\|
[z16.a0g0b2i875.0.zip](ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/Ip-BETA/z16.a0g0b2i875.0.zip)
\|\|
[z16.a0g0b3i875.0.zip](ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/Ip-BETA/z16.a0g0b3i875.0.zip)
\|\|
[z16.a0g0b4i875.1.zip](ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/Ip-BETA/z16.a0g0b4i875.1.zip)
\|\| XX \|\| XX \|\| \|\| 131 \|\|
[z16.a0g0b0i131.3.zip](ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/Ip-BETA/z16.a0g0b0i131.3.zip)
\|\|
[z16.a0g0b1i131.a.zip](ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/Ip-BETA/z16.a0g0b1i131.a.zip)
\|\|
[z16.a0g0b2i131.2.zip](ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/Ip-BETA/z16.a0g0b2i131.2.zip)
\|\|
[z16.a0g0b3i131.1.zip](ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/Ip-BETA/z16.a0g0b3i131.1.zip)
\|\|
[z16.a0g0b4i131.4.zip](ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/Ip-BETA/z16.a0g0b4i131.4.zip)
\|\| XX \|\| XX \|\| \|\| 174 \|\|
[z16.a0g0b0.5.zip](ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/Ip-BETA/z16.a0g0b0.5.zip)
\|\|
[z16.a0g0b1.4.zip](ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/Ip-BETA/z16.a0g0b1.4.zip)
\|\|
[z16.a0g0b2.2.zip](ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/Ip-BETA/z16.a0g0b2.2.zip)
\|\|
[z16.a0g0b3.3.zip](ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/Ip-BETA/z16.a0g0b3.3.zip)
\|\|
[z16.s3.1.zip](ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/Ip-BETA/z16.s3.1.zip)
\|\|
[z16.a0g0b5.5.zip](ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/Ip-BETA/z16.a0g0b5.5.zip)
\|\| XX \|\| \|\| 200 \|\| XX \|\| XX \|\| XX \|\| XX \|\| XX \|\| XX
\|\|
[z16.a0g0b6.9b.zip](ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/Ip-BETA/z16.a0g0b6.9b.zip)
\|\|

### Core Current Profile Variation (Table VII)

\|\| Alpha \|\| File Location \|\| \|\| 0.0 \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/CURRENT.PROFILE.SCAN/z16.a0g0b3.3.zip>
\|\| \|\| 0.1 \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/CURRENT.PROFILE.SCAN/z16.a1g0b3.0.zip>
\|\| \|\| 0.2 \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/CURRENT.PROFILE.SCAN/z16.a2g0b3.0.zip>
\|\| \|\| 0.3 \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/CURRENT.PROFILE.SCAN/z16.a3g0b3.0.zip>
\|\| \|\| 0.4 \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/CURRENT.PROFILE.SCAN/z16.a4g0b3.1.zip>
\|\| \|\| 0.5 \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/CURRENT.PROFILE.SCAN/z16.a5g0b3.3.zip>
\|\|

### Edge Current Profile Variation

\|\| Delta \|\| File Location \|\| \|\| 0.0 \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/CURRENT.PROFILE.SCAN/z16.a0g0b5.5.zip>
\|\| \|\| 0.2 \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/CURRENT.PROFILE.SCAN/z16.a0g0b5.5.min+Jedge0.2.zip>
\|\| \|\| 0.4 \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/CURRENT.PROFILE.SCAN/z16.a0g0b5.5.min+Jedge0.4.zip>
\|\|

### Core Pressure Profile Variation (Table X)

\|\| Gamma \|\| File Location \|\| \|\| 0.0 \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/PRESSURE.PROFILE.SCAN/z16.a0g0b3.3.zip>
\|\| \|\| 0.2 \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/PRESSURE.PROFILE.SCAN/z16.a0g2b3.a.zip>
\|\| \|\| 0.4 \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/PRESSURE.PROFILE.SCAN/z16.a0g4b3.a.zip>
\|\| \|\| 0.6 \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/PRESSURE.PROFILE.SCAN/z16.a0g6b3.b.zip>
\|\| \|\| 0.8 \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/PRESSURE.PROFILE.SCAN/z16.a0g8b3.c.zip>
\|\|

### Pressure Pedestal

<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/PRESSURE.PROFILE.SCAN/z16.a0b3ppedestal.3c.zip>

### Variation of Iota at Fixed Shear (Table XII)

\|\| Delta Iota \|\| File Location \|\| \|\| Ref \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/S3.IOTA.SCAN/z16.s3ioscan.1.zip>
\|\| \|\| -0.1 \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/S3.IOTA.SCAN/z16.s3ioscan.3.zip>
\|\| \|\| -0.2 \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/S3.IOTA.SCAN/z16.s3ioscan.5c.zip>
\|\| \|\| +0.1 \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/S3.IOTA.SCAN/z16.s3ioscan.2a.zip>
\|\|

### Variation of Shear at Fixed Core Iota

\|\| Edge Iota \|\| File Location \|\| \|\| 0.58 \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/S3.IOTA.SCAN/0330h.zip>
\|\| \|\| 0.64 (ref) \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/S3.IOTA.SCAN/z16.s3ioscan.1.zip>
\|\| \|\| 0.74 \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/S3.IOTA.SCAN/z16.s3ioscan.7a.zip>
\|\| \|\| 0.84 \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/S3.IOTA.SCAN/0319b.zip>
\|\|

### Shape Stabilization Experiments (Table XIV)

\|\| Experiment \|\| File Location \|\| \|\| A \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/3DSHAPE-STABILIZATION/0420a.zip>
\|\| \|\| B \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/3DSHAPE-STABILIZATION/0413b.zip>
\|\| \|\| C \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/3DSHAPE-STABILIZATION/0412e.zip>
\|\|

### Neoclassical Transport Variation Study (Table XV)

\|\| Optimization Level \|\| File Location \|\| \|\| Fully Degraded \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/QA.SCAN/0325qa.2.zip>
\|\| \|\| Partially Degraded \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/QA.SCAN/0325qa.3.zip>
\|\| \|\| Reference Optimized \|\|
<ftp://ftp.pppl.gov/pub/pomphrey/M45_1231k10j2/QA.SCAN/z16.a0g0b2i875.0.zip>
\|\|
