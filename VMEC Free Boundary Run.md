Tutorial: VMEC Free Boundary Run
================================

This tutorial will walk the user through running VMEC with a free
boundary condition. For this example the National Compact Stellarator
Experiment (NCSX) configuration will be used. This machine is
stellarator symmetric with a periodicity of three.

------------------------------------------------------------------------

## Create the 'mgrid' file 

First one must use the
[MAKEGRID](MAKEGRID) (XGRID) routine to process the coils definition
file and produce an 'mgrid' file. For this example we will use the
'coils.c09r00' file which has been developed for the NCSX device
([coils.c09r00](examples/coils.c09r00)). We create a text file
'input_xgrid.dat' which contains:

    c09r00
    S
    y
    0.436
    2.436
    -1.0
    1.0
    36
    201
    201

Note that line two should be omitted for older (pre 8.0 versions of
the VMEC suite of codes). > The xgrid code may then be executed with
this input file (this may take some time depending on your computing
power).

    > ~/bin/xgrid < input_xgrid.dat >& log_xgrid.c09r00 &

or type it interactively

     Enter extension of "coils" file     :  
     Scale (S) bfield to unit current/turn OR use raw (R) currents from coils file:  
     Assume stellarator symmetry (Y/N)?  :  
     Enter rmin (min radial grid dimension)  :  
     Enter rmax (max radial grid dimension)  :  
     Enter zmin (min vertical grid dimension):  
     Enter zmax (max vertical grid dimension):  
     Enter number of toroidal planes/period  :  
     Enter number of r (radial) mesh points  :  
     Enter number of z mesh points  : 

The code will begin to run and part of the screen output is:

     Stellarator symmetry IS assumed
     rmin =  0.436  rmax =  2.4359999999999999
     zmin =  -1.  zmax =  1.
     kp =  36  ir =  201  jz =  201
     
     Input  file: coils.c09r00
     Mgrid  file: mgrid_c09r00
     Extcur file: extcur.c09r00
     COIL GROUP          :  ModA
     TOTAL COILS IN GROUP:   6 TOTAL FILAMENTS:   2400
     K =    1 (OUT OF   36)
     K =    2
     K =    3
     K =    4
     K =    5
     K =    6
     K =    7
     K =    8
     K =    9
     K =   10
     K =   11
     K =   12
     K =   13
     K =   14
     K =   15
     K =   16
     K =   17
     K =   18
     K =   19
     COIL GROUP          :  ModB
     TOTAL COILS IN GROUP:   6 TOTAL FILAMENTS:   2400
     K =    1 (OUT OF   36)
     K =    2
     
    .
    .
    .
     COIL GROUP          :  TF
     TOTAL COILS IN GROUP:  18 TOTAL FILAMENTS:   2196
     K =    1 (OUT OF   36)
     K =    2
     K =    3
     K =    4
     K =    5
     K =    6
     K =    7
     K =    8
     K =    9
     K =   10
     K =   11
     K =   12
     K =   13
     K =   14
     K =   15
     K =   16
     K =   17
     K =   18
     K =   19
     TIME IN PARSER =    0.397 SECONDS
     TIME IN BFIELD =  407.148 SECONDS
     THE BFIELDS HAVE BEEN STORED IN THE MGRID FILE IN SCALED MODE. THE EXTERNAL
     CURRENTS CORRESPONDING TO THOSE IN THE COILS-DOT FILE
     ARE GIVEN IN THE EXTCUR ARRAY IN THE FILE extcur.c09r00. 
     THEY SHOULD BE ENTERED INTO THE VMEC INPUT (INDATA) FILE.

This choice of parameters gives us [cm] scale resolution and
provides us with enough toroidal resolution to run with 8 toroidal modes
per field period. At the end of the run you should have produced the
following files: `mgrid_c09r00.nc` (or `mgrid.c09r00`), `extcur.c09r00`, and
`log_xgrid.c09r00`.  

## Edit the input namelist text file.

The input namelist
([input.ncsx_c09r00_free](examples/input.ncsx_c09r00_free))
controls the execution of the VMEC code. The suffix of the input file
will be appended to each of the output files as we will see after
execution. The Fourier coefficient in this file have been generated
through an optimization routine. In general, more simple initial
conditions will suffice for the axis position and outer most flux
surface. The name of the mgrid file must now be specified, along with
setting `LFREEB` to true. The `NZETA` variable must match the value you
selected for mgrid creation (see 7th line in the `input_xgrid.dat`
file). You'll notice that for the free boundary run the user must also
specify the EXTCUR array. This array specifies the current running
through each coil group. There should be one entry per coil group. The
`extcur` file contains suggestions for these values based on the
values it read from the `coils` file. These values can be copied and
pasted directly into the VMEC 'input' file. Note that the traditional
polynomial form of the current profile (NCURR=1) and pressure profile
are being used. Note that if you have a binary 'mgrid' file then
you'll need to modify the MGRID_FILE variable to match the proper
name.

    &INDATA
    !----- Runtime Parameters -----
      DELT =   9.00E-01
      NITER = 5000
      NSTEP =  200
      TCON0 =   2.00E+00
      NS_ARRAY =   9   29   49   99
      FTOL_ARRAY = 1.000000E-06  1.000000E-08  1.000000E-10  1.000000E-12
    !----- Grid Parameters -----
      LASYM = F
      NFP =    3
      MPOL =   11
      NTOR =    6
      NZETA  =   36
      PHIEDGE =   4.97070205837336E-01
    !----- Free Boundary Parameters -----
      LFREEB = T
      MGRID_FILE = 'mgrid_c09r00.nc'
      EXTCUR =   6.52271941985300E+05  6.51868569367400E+05  5.37743588647300E+05
      2.50000000000000E-07  2.50000000000000E-07  2.80949750000000E+04
      -5.48049500000000E+04  3.01228950000000E+04  9.42409100000000E+04
      4.55138737653200E+04
      NVACSKIP =    6
    !----- Pressure Parameters -----
      GAMMA =   0.000000E+00
      BLOAT =   1.000000E+00
      SPRES_PED =   1.00000000000000E+00
      AM =   6.85517649352426E+04 -5.12027745123057E+03 -3.61510451745464E+04 -4.74263014113066E+05
      1.78878195473870E+06 -3.21513828868170E+06  2.69041023837233E+06 -8.17049854168367E+05
      0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
    !----- Current/Iota Parameters -----
      CURTOR =  -1.78606250000000E+05
      NCURR =    1
      AI =   0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
      0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
      0.00000000000000E+00  0.00000000000000E+00  0.00000000000000E+00
      AC =   8.18395699999999E+03  1.43603560000000E+06 -1.07407140000000E+07  7.44389200000000E+07
     -3.22215650000000E+08  8.81050800000000E+08 -1.49389660000000E+09  1.52746800000000E+09
     -8.67901590000000E+08  2.10351200000000E+08  0.00000000000000E+00
    !----- Axis Parameters -----
      RAXIS =   1.49569454253276E+00  1.05806400912320E-01  7.21255454715878E-03 -3.87402652289249E-04
     -2.02425864534069E-04 -1.62602353744308E-04 -8.89569831063077E-06
      ZAXIS =   0.00000000000000E+00 -5.19492027001782E-02 -3.18814224021375E-03  2.26199929262002E-04
      1.28803681387330E-04  1.11266150452637E-06  1.13732703961869E-05
    !----- Boundary Parameters -----
      RBC(0,0) =   1.40941668895656E+00     ZBS(0,0) =   0.00000000000000E+00
      RBC(1,0) =   2.79226697269945E-02     ZBS(1,0) =  -1.92433268059631E-02
      RBC(2,0) =  -1.54739398509667E-03     ZBS(2,0) =   1.11459511078088E-02
      RBC(3,0) =   2.90733840104882E-03     ZBS(3,0) =  -3.97869471888770E-03
      RBC(4,0) =  -8.91322016448873E-04     ZBS(4,0) =   1.34394804673514E-03
      RBC(5,0) =  -7.81997770407839E-05     ZBS(5,0) =  -1.57143910159387E-04
      RBC(6,0) =   1.06129711928351E-04     ZBS(6,0) =   9.58024291307491E-05
      RBC(-6,1) =  2.48228899767757E-05     ZBS(-6,1) = -2.28386224209054E-05
      RBC(-5,1) =  8.23567667077671E-05     ZBS(-5,1) =  3.30176003890210E-04
      RBC(-4,1) = -7.20220898033597E-04     ZBS(-4,1) =  1.28038328362904E-04
      RBC(-3,1) =  2.76250777733235E-03     ZBS(-3,1) =  3.43199911886317E-04
      RBC(-2,1) = -1.24883373588382E-02     ZBS(-2,1) =  6.12174680232785E-04
      RBC(-1,1) =  1.52272804511910E-02     ZBS(-1,1) = -2.70066914159594E-02
      RBC(0,1) =   2.89195233044040E-01     ZBS(0,1) =   4.50462554508443E-01
      RBC(1,1) =  -1.17988850341728E-01     ZBS(1,1) =   1.93490230971634E-01
      RBC(2,1) =  -3.84923299492945E-03     ZBS(2,1) =   5.72865331625290E-03
      RBC(3,1) =  -1.44452305429529E-03     ZBS(3,1) =   2.19788951889214E-03
      RBC(4,1) =  -2.11622985211109E-04     ZBS(4,1) =   1.31883972780290E-03
      RBC(5,1) =   1.79091719677667E-04     ZBS(5,1) =  -5.63363462408534E-04
      RBC(6,1) =   1.31982402741742E-04     ZBS(6,1) =  -9.31801467009349E-05
      RBC(-6,2) = -2.40882614870476E-05     ZBS(-6,2) = -3.95416405717970E-05
      RBC(-5,2) = -4.92449386382591E-05     ZBS(-5,2) = -3.25048356502217E-06
      RBC(-4,2) =  1.50530476034115E-04     ZBS(-4,2) =  4.61522421935086E-05
      RBC(-3,2) = -1.23084235126550E-03     ZBS(-3,2) = -3.40868203306282E-04
      RBC(-2,2) =  2.01350576071929E-04     ZBS(-2,2) = -4.19781517712033E-03
      RBC(-1,2) =  2.36777003797179E-03     ZBS(-1,2) =  1.98753868216412E-02
      RBC(0,2) =   5.73443941583452E-02     ZBS(0,2) =   4.81527027892127E-03
      RBC(1,2) =   6.89385874058265E-02     ZBS(1,2) =  -9.28353553039424E-03
      RBC(2,2) =   4.71996849673782E-02     ZBS(2,2) =  -2.04292782322197E-02
      RBC(3,2) =  -5.50889052720066E-04     ZBS(3,2) =   8.81593501270446E-04
      RBC(4,2) =   4.24491391207156E-04     ZBS(4,2) =  -6.08871281835245E-04
      RBC(5,2) =  -2.07538883155595E-04     ZBS(5,2) =  -3.88708113241096E-04
      RBC(6,2) =  -1.62304038006678E-04     ZBS(6,2) =   1.72340342752605E-04
      RBC(-6,3) = -1.01105699684233E-04     ZBS(-6,3) = -6.16215454248342E-05
      RBC(-5,3) =  5.15925605980565E-05     ZBS(-5,3) =  1.23419431936950E-04
      RBC(-4,3) = -3.79290487874111E-05     ZBS(-4,3) =  3.98637008165582E-06
      RBC(-3,3) = -2.96154201246223E-04     ZBS(-3,3) = -7.01248486620889E-04
      RBC(-2,3) =  1.27628943631957E-03     ZBS(-2,3) =  3.19332333533202E-03
      RBC(-1,3) =  3.12803506573940E-03     ZBS(-1,3) = -8.24657727838880E-03
      RBC(0,3) =  -1.34574092972690E-02     ZBS(0,3) =   5.05936199755365E-03
      RBC(1,3) =  -8.02339287294677E-03     ZBS(1,3) =  -3.90421394288867E-03
      RBC(2,3) =  -1.68510947837154E-02     ZBS(2,3) =   3.75441853342170E-03
      RBC(3,3) =  -8.00581733372124E-03     ZBS(3,3) =   6.00542774606014E-03
      RBC(4,3) =   1.80667899211621E-03     ZBS(4,3) =  -4.16787432635077E-04
      RBC(5,3) =   3.10773970094350E-05     ZBS(5,3) =   5.44335921432213E-05
      RBC(6,3) =   8.32496816115997E-05     ZBS(6,3) =  -4.15830451164888E-05
      RBC(-6,4) = -1.19874891436340E-05     ZBS(-6,4) =  1.56845408711308E-05
      RBC(-5,4) =  1.22793444338155E-04     ZBS(-5,4) = -3.97576733690054E-05
      RBC(-4,4) = -1.30945484439682E-04     ZBS(-4,4) = -7.22429623460448E-05
      RBC(-3,4) = -1.21368603604647E-04     ZBS(-3,4) =  3.52928331257216E-04
      RBC(-2,4) =  1.00352526472782E-03     ZBS(-2,4) = -1.23710282249961E-04
      RBC(-1,4) = -1.73680844498789E-03     ZBS(-1,4) = -1.50689928334813E-03
      RBC(0,4) =   1.80149787198970E-03     ZBS(0,4) =   1.56109492686192E-03
      RBC(1,4) =   3.82771889154294E-03     ZBS(1,4) =   3.80910842862487E-03
      RBC(2,4) =   5.43835034437129E-03     ZBS(2,4) =   2.06275075117804E-03
      RBC(3,4) =   8.39729828422411E-04     ZBS(3,4) =  -1.54779126563731E-03
      RBC(4,4) =   6.74263596810560E-04     ZBS(4,4) =  -1.33149943553452E-03
      RBC(5,4) =  -6.98647584180715E-04     ZBS(5,4) =   3.81307095116973E-04
      RBC(6,4) =   8.77670652920776E-05     ZBS(6,4) =  -1.40433963574141E-05
      RBC(-6,5) =  6.78635213884316E-06     ZBS(-6,5) = -1.22283666932084E-05
      RBC(-5,5) =  3.87846546342867E-05     ZBS(-5,5) =  4.64829761643373E-05
      RBC(-4,5) = -3.78300368387435E-05     ZBS(-4,5) = -7.03801581329045E-05
      RBC(-3,5) = -1.21743926248229E-05     ZBS(-3,5) =  1.85735151533626E-04
      RBC(-2,5) = -2.68229697014545E-04     ZBS(-2,5) = -9.33216243296025E-04
      RBC(-1,5) =  1.19567316567517E-03     ZBS(-1,5) =  2.12648562837673E-03
      RBC(0,5) =  -7.12579133390599E-04     ZBS(0,5) =  -1.97890515574565E-03
      RBC(1,5) =   8.81127157923892E-04     ZBS(1,5) =   2.71321673191593E-03
      RBC(2,5) =   9.67210453659238E-04     ZBS(2,5) =   8.74618447862515E-04
      RBC(3,5) =   2.11794179698155E-04     ZBS(3,5) =   8.43817701627930E-04
      RBC(4,5) =   1.29403911922840E-03     ZBS(4,5) =   6.51808476607835E-04
      RBC(5,5) =  -1.30477683585083E-04     ZBS(5,5) =   1.01349326961770E-04
      RBC(6,5) =   1.86680624010370E-04     ZBS(6,5) =  -2.13838628730300E-04
      RBC(-6,6) = -4.08213549686361E-05     ZBS(-6,6) = -7.53394429655583E-06
      RBC(-5,6) =  7.11305157811999E-05     ZBS(-5,6) =  2.54876062250879E-05
      RBC(-4,6) = -1.33727065581923E-04     ZBS(-4,6) = -1.70180862196520E-05
      RBC(-3,6) =  1.65191943182183E-06     ZBS(-3,6) = -1.31350577800873E-04
      RBC(-2,6) =  2.19460449719541E-04     ZBS(-2,6) =  4.38914760402648E-04
      RBC(-1,6) =  4.68618562605432E-04     ZBS(-1,6) = -4.44537659614533E-04
      RBC(0,6) =  -8.51896573200937E-04     ZBS(0,6) =   7.36122964253313E-04
      RBC(1,6) =  -5.26623264534578E-05     ZBS(1,6) =  -1.12352425125337E-03
      RBC(2,6) =  -1.31954654361710E-04     ZBS(2,6) =  -2.22905186553194E-03
      RBC(3,6) =  -8.91482312658694E-04     ZBS(3,6) =  -2.11193996461398E-03
      RBC(4,6) =  -3.89733094884781E-04     ZBS(4,6) =  -3.44184359663702E-04
      RBC(5,6) =  -2.74329775462215E-04     ZBS(5,6) =  -5.06914660659672E-05
      RBC(6,6) =   2.47385092660320E-04     ZBS(6,6) =   3.74971583066409E-05
      RBC(-6,7) =  9.61516193308531E-06     ZBS(-6,7) = -3.66121037894761E-06
      RBC(-5,7) = -2.51122684780459E-05     ZBS(-5,7) =  3.72828134065079E-05
      RBC(-4,7) =  4.44568599556351E-05     ZBS(-4,7) = -8.74488353626824E-05
      RBC(-3,7) = -1.42433799354752E-04     ZBS(-3,7) =  1.48694485468843E-04
      RBC(-2,7) =  4.85802385952487E-04     ZBS(-2,7) = -2.27519962800893E-04
      RBC(-1,7) = -9.00652688032426E-04     ZBS(-1,7) =  4.16601324903870E-04
      RBC(0,7) =   9.59457670863182E-04     ZBS(0,7) =  -3.25818663499641E-04
      RBC(1,7) =  -3.37159659594826E-04     ZBS(1,7) =  -2.34240245561361E-04
      RBC(2,7) =  -4.64969900861713E-04     ZBS(2,7) =   4.87821281121050E-04
      RBC(3,7) =  -4.09185322970312E-04     ZBS(3,7) =   8.50140634573578E-04
      RBC(4,7) =   5.32088748759921E-05     ZBS(4,7) =   5.93528572346752E-04
      RBC(5,7) =  -3.21692982976907E-04     ZBS(5,7) =  -2.54775193277671E-04
      RBC(6,7) =  -4.82403633897412E-05     ZBS(6,7) =   1.41947169759239E-05
      RBC(-6,8) = -2.23522770283961E-05     ZBS(-6,8) = -4.00911971000495E-06
      RBC(-5,8) =  3.95696912099304E-05     ZBS(-5,8) =  1.34684147523625E-05
      RBC(-4,8) = -6.50775924544567E-05     ZBS(-4,8) = -2.94168940555405E-05
      RBC(-3,8) =  1.71610112932980E-04     ZBS(-3,8) =  2.17875987311858E-05
      RBC(-2,8) = -3.45412623614909E-04     ZBS(-2,8) = -7.26482153663716E-05
      RBC(-1,8) =  5.61089095467387E-04     ZBS(-1,8) =  2.51145295676537E-04
      RBC(0,8) =  -5.84359101746051E-04     ZBS(0,8) =  -5.42465826224607E-04
      RBC(1,8) =  -6.16860761080513E-05     ZBS(1,8) =   3.93697603313273E-04
      RBC(2,8) =   5.99275780897287E-04     ZBS(2,8) =   3.30798770955874E-04
      RBC(3,8) =   5.68520162541870E-04     ZBS(3,8) =   5.47788467933391E-04
      RBC(4,8) =   4.47404034542356E-04     ZBS(4,8) =   2.43547539548605E-04
      RBC(5,8) =   2.76704814165950E-04     ZBS(5,8) =   9.15194583315619E-05
      RBC(6,8) =   2.97621090888441E-04     ZBS(6,8) =   1.65605427353701E-04
      RBC(-6,9) = -3.78145897931544E-06     ZBS(-6,9) = -1.85759364750771E-06
      RBC(-5,9) = -1.57985677623482E-05     ZBS(-5,9) =  1.06970045147331E-05
      RBC(-4,9) =  7.91641381274532E-05     ZBS(-4,9) = -1.16252171939772E-05
      RBC(-3,9) = -1.97587428419659E-04     ZBS(-3,9) = -3.08457797690412E-05
      RBC(-2,9) =  3.95855751452672E-04     ZBS(-2,9) =  2.03418980231168E-04
      RBC(-1,9) = -5.41153103221438E-04     ZBS(-1,9) = -3.99552958537408E-04
      RBC(0,9) =   4.98714381092541E-04     ZBS(0,9) =   4.32916759100965E-04
      RBC(1,9) =  -8.06048953531492E-05     ZBS(1,9) =  -1.84722027458208E-04
      RBC(2,9) =  -8.67990109801738E-05     ZBS(2,9) =   2.52568631885491E-04
      RBC(3,9) =   4.35340840113358E-04     ZBS(3,9) =  -3.50159782847442E-04
      RBC(4,9) =   2.33585243788111E-04     ZBS(4,9) =  -7.06133299107118E-04
      RBC(5,9) =  -7.69581174305243E-06     ZBS(5,9) =  -3.79072907220561E-04
      RBC(6,9) =  -2.85256407945938E-05     ZBS(6,9) =  -4.49599333610498E-05
      RBC(-6,10) =  1.20206720198758E-05     ZBS(-6,10) =  4.73580005255806E-06
      RBC(-5,10) =  7.02670536357846E-06     ZBS(-5,10) = -6.99664911015022E-06
      RBC(-4,10) = -2.76926398374910E-05     ZBS(-4,10) = -9.18014408856618E-06
      RBC(-3,10) =  5.20223745639364E-05     ZBS(-3,10) =  6.80574180383753E-05
      RBC(-2,10) = -7.88310431749746E-05     ZBS(-2,10) = -1.06370673487973E-04
      RBC(-1,10) =  1.21755712542490E-05     ZBS(-1,10) =  1.22161894513591E-04
      RBC(0,10) =  3.22193519645521E-05     ZBS(0,10) = -6.04052049877600E-05
      RBC(1,10) = -1.08453911913102E-06     ZBS(1,10) =  8.60890353665843E-05
      RBC(2,10) =  1.04051545504927E-04     ZBS(2,10) = -2.17661420286656E-04
      RBC(3,10) = -5.21965328013036E-04     ZBS(3,10) = -2.67111216700977E-04
      RBC(4,10) = -4.95991087393098E-04     ZBS(4,10) =  2.43875640076056E-05
      RBC(5,10) = -1.94520415280627E-04     ZBS(5,10) =  1.55759001593971E-04
      RBC(6,10) = -6.94143617569942E-05     ZBS(6,10) =  4.40565098025554E-05
    /
    &END

## Execute the code.

To simplify execution of the code,
the VMEC compilation scripts create a directory called `bin` in your
home () directory. Symbolic links are then placed there pointing to each
of the compiled codes in their respective 'Vrelease' subdirectories.
In practice, the screen output from VMEC should be redirected to a log
file and put in the background (>& log.ncsx_c09r00_free &). This is
done by passing the suffix of the input file to the VMEC code through
the command line.

    >~/bin/xvmec2000 ncsx_c09r00_free
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      VMEC OUTPUT FILES ALREADY EXIST: OVERWRITING THEM ...
      SEQ =    1 TIME SLICE  0.0000E+00
      PROCESSING INPUT.ncsx_c09r00_free
      THIS IS VMEC2000, A 3D EQUILIBRIUM CODE, VERSION 8.47
      Lambda: Full Radial Mesh. L-Force: hybrid full/half.
      
      COMPUTER: computer.domain.net   OS: Linux   RELEASE: 2.6.18-194.17.4.el5  DATE = Sep 06,2011  TIME = 13:13:44
      Opening vacuum field file: mgrid_c09r00.nc
      Time to read MGRID file:   3.10E+00 s
      NS =    9 NO. FOURIER MODES =  137 FTOLV =  1.000E-06 NITER =   5000
     ITER    FSQR      FSQZ      FSQL     RAX(v=0)      WMHD      DEL-BSQ
        1  5.24E+01  6.99E+00  1.63E-01  1.608E+00  3.7781E+00  1.000E+00
      In VACUUM, np =  3  mf = 12  nf =  6 nu = 28  nv =   36
      2*pi * a * -BPOL(vac) =  -1.89E-01 TOROIDAL CURRENT =  -1.89E-01
      R * BTOR(vac) =   2.37E+00 R-BTOR =   2.38E+00
      VACUUM PRESSURE TURNED ON AT   82 ITERATIONS
      165  9.55E-07  2.25E-07  4.60E-07  1.609E+00  3.5227E+00  4.634E-03
      NS =   29 NO. FOURIER MODES =  137 FTOLV =  1.000E-08 NITER =   5000
     ITER    FSQR      FSQZ      FSQL     RAX(v=0)      WMHD      DEL-BSQ
        1  6.96E-02  3.27E-02  3.97E-04  1.609E+00  3.5227E+00  4.634E-03
      200  6.23E-08  1.59E-08  8.03E-09  1.609E+00  3.5199E+00  2.732E-03
      287  9.99E-09  2.61E-09  2.17E-09  1.609E+00  3.5199E+00  2.818E-03
      NS =   49 NO. FOURIER MODES =  137 FTOLV =  1.000E-10 NITER =   5000
     ITER    FSQR      FSQZ      FSQL     RAX(v=0)      WMHD      DEL-BSQ
        1  3.14E+00  1.57E+00  4.90E-06  1.609E+00  3.5198E+00  2.818E-03
      200  2.19E-08  8.35E-09  1.21E-09  1.608E+00  3.5197E+00  2.585E-03
      400  1.28E-09  3.24E-10  2.86E-10  1.608E+00  3.5197E+00  2.756E-03
      600  3.28E-10  8.59E-11  5.96E-11  1.608E+00  3.5197E+00  2.848E-03
      782  9.95E-11  2.66E-11  1.57E-11  1.608E+00  3.5197E+00  2.873E-03
      NS =   99 NO. FOURIER MODES =  137 FTOLV =  1.000E-12 NITER =   5000
     ITER    FSQR      FSQZ      FSQL     RAX(v=0)      WMHD      DEL-BSQ
        1  9.60E+00  4.73E+00  2.02E-06  1.608E+00  3.5197E+00  2.873E-03
      200  3.97E-08  2.25E-08  9.27E-11  1.608E+00  3.5196E+00  2.605E-03
      400  1.13E-09  1.52E-10  4.09E-11  1.608E+00  3.5196E+00  2.611E-03
      600  3.86E-10  5.73E-11  2.43E-11  1.608E+00  3.5196E+00  2.620E-03
      800  1.49E-10  2.82E-11  1.36E-11  1.608E+00  3.5196E+00  2.627E-03
     1000  5.74E-11  1.33E-11  6.58E-12  1.608E+00  3.5196E+00  2.632E-03
     1200  2.00E-11  4.94E-12  2.76E-12  1.608E+00  3.5196E+00  2.634E-03
     1400  7.79E-12  1.87E-12  1.13E-12  1.608E+00  3.5196E+00  2.635E-03
     1600  3.57E-12  8.99E-13  4.91E-13  1.608E+00  3.5196E+00  2.634E-03
     1800  1.39E-12  3.28E-13  2.27E-13  1.608E+00  3.5196E+00  2.634E-03
     1881  9.93E-13  2.53E-13  1.68E-13  1.608E+00  3.5196E+00  2.634E-03
     EXECUTION TERMINATED NORMALLY
     FILE : ncsx_c09r00_free
     NUMBER OF JACOBIAN RESETS =    3
        TOTAL COMPUTATIONAL TIME             501.04 SECONDS
        TIME TO READ IN DATA                   3.13 SECONDS
        TIME TO WRITE DATA TO WOUT             0.30 SECONDS
        TIME IN EQFORCE                        7.22 SECONDS
        TIME IN VACUUM LOOP                   82.27 SECONDS
        TIME IN FOURIER TRANSFORM            117.16 SECONDS
        TIME IN INVERSE FOURIER XFORM         80.03 SECONDS
        TIME IN FORCES                        95.99 SECONDS
        TIME IN BCOVAR                        77.08 SECONDS
        TIME IN RESIDUE                        2.29 SECONDS
        TIME (REMAINDER) IN FUNCT3D           34.30 SECONDS

## Examine the output.
 
For this example four files were
created (`jxbout.ncsx_c09r00_free`, `mercier.ncsx_c09r00_free`,
`threed1.ncsx_c09r00_free`, and `wout.ncsx_c09r00_free`). As was
mentioned before, each file had the suffix of the input file appended to
it's name. This allows multiple runs to be stored in the same directory
for comparison. The `jxbout` file contains values for various
quantities on a grid throughout the simulation domain. The `mercier`
file contains radial profiles (radial index in VMEC is denoted by the
variable `s`) of various quantities. The `threed1` file can be
considered an expanded log file where various quantities are calculated
which were not output to the screen. This file is fairly self
explanatory. The `wout` file is the data file for the run. It contains
the Fourier Coefficients for the magnetic field along with various
quantities. A few packages exist to visualize this data and the user is
encourage to use these as templates for their own visualization
routines.
