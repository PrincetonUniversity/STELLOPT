#!/usr/bin/env python3
import sys, os
sys.path.insert(0, '../../../pySTEL/')
import numpy as np                    #For Arrays
from math import pi
from libstell.stellopt import read_stellopt

try:
	qtCreatorPath=os.environ["STELLOPT_PATH"]
except KeyError:
	print("Please set environment variable STELLOPT_PATH")
	sys.exit(1)

failtol = 5.0
filename='stellopt.BASIC'
data=read_stellopt(filename)
if not data:
    print('ERROR Opening File: '+filename)
    sys.exit(0)
else:
    print('EXTENSION: '+filename)
print('==== Scalars ====')
varlist={}
varlist['PHIEDGE_equil']=0.514386
varlist['CURTOR_equil']=-1.742957864085E+005
varlist['RBTOR_equil']=2.312873858830
varlist['R0_equil']=1.577564814144
varlist['VOLUME_equil']=2.9787172145375
varlist['BETA_equil']=4.247495078898E-002
varlist['WP_equil']=1.925493826145E+005
varlist['ASPECT_equil']=4.365254725961
lfail = 0;
for temp in varlist:
    act = varlist[temp]
    cal = data[temp]
    perct = 100*(abs(act-cal)/act)
    print('  '+temp+': '+str(cal)+'   '+str(act)+'   '+str(int(perct))+'%')
    if perct > failtol:
        lfail = 1
print('==== Vectors ====')
varlist={}

varlist['NELINE_equil']=np.array([3.953450576613e+19, 4.3636068324e+19, 4.991048247502e+19, 4.971516997227e+19, 4.951985746951e+19, 4.989420643313e+19, 4.960123767899e+19, 4.989420643313e+19, 4.950358142762e+19, 4.969889393037e+19, 4.979655018175e+19, 4.989420643313e+19, 4.989420643313e+19, 4.979655018175e+19, 4.961751372089e+19, 4.932454496676e+19, 4.969889393037e+19, 4.991048247503e+19, 4.9129232464e+19, 4.991048247502e+19, 4.950358142762e+19, 4.969889393037e+19, 4.940592517624e+19, 4.951985746951e+19, 4.991048247502e+19, 4.922688871538e+19, 4.979655018175e+19, 0.0, 0.0, 0.0])
varlist['FARADAY_equil']=np.array([4.894751765575e+18, 6.372097940318e+18, 8.439039155694e+18, 8.164499074463e+18, 7.907239879865e+18, 7.658738541811e+18, 7.417446825527e+18, 7.185683021612e+18, 6.961912643751e+18, 6.743210177623e+18, 6.523837515289e+18, 6.30721842952e+18, 6.090569474666e+18, 5.890734466108e+18, 5.698670754699e+18, 5.509665279525e+18, 5.332682574496e+18, 5.163807294034e+18, 5.012346549374e+18, 4.868832246336e+18, 4.727762963751e+18, 4.572504179913e+18, 4.390247914634e+18, 4.138238855357e+18, 3.744489625987e+18, 3.026773488722e+18, 1.512169572346e+18, 0.0, 0.0, 0.0])
varlist['PRESS_equil']=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 261.1067539067, 7567.128185418, 20047.14655471, 32312.31159263, 42332.23222524, 50046.13700725, 55917.42539314, 60408.36983676, 63860.68080978, 66505.23749345, 68502.24705172, 69975.56073848, 71031.46403811, 71765.1809579, 72260.29198991, 72586.37376406, 72796.91554, 72930.46126609, 73012.03896046, 73057.41843957, 73076.67616391, 73076.21349688, 73061.55757927, 73037.69093725, 73009.08243876, 72979.47978793, 72952.73759202, 72932.50276448, 72921.60815253, 72921.60819293, 72932.50288191, 72952.73777575, 72979.48002199, 73009.08270461, 73037.69122059, 73061.55786661, 73076.21377578, 73076.67642487, 73057.41867484, 73012.03916128, 72930.4614214, 72796.9156205, 72586.37371961, 72260.29174039, 71765.18038716, 71031.46299259, 69975.55902886, 68502.24445283, 66505.23374114, 63860.67557666, 60408.36267824, 55917.41565643, 50046.12373057, 42332.21415052, 32312.28762336, 20047.11769959, 7567.102407716, 261.1007023661, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
varlist['NE_equil']=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
varlist['TE_equil']=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 15.72535141048, 471.9346609297, 1251.161798764, 2016.861833142, 2642.328553287, 3123.801034882, 3490.246584935, 3770.537216928, 3986.010840727, 4151.075868028, 4275.729194636, 4367.694630125, 4433.605235154, 4479.400532239, 4510.29898452, 4530.640380848, 4543.776921854, 4552.099612378, 4557.190921618, 4560.035324749, 4561.236074888, 4561.207243668, 4560.29358245, 4558.80212399, 4557.00527658, 4555.153708361, 4553.487281259, 4552.226843861, 4551.547617573, 4551.547620092, 4552.226851179, 4553.487292701, 4555.153722962, 4557.005293266, 4558.802141736, 4560.293600374, 4561.207261049, 4561.23609115, 4560.035339433, 4557.190934234, 4552.099622058, 4543.776926867, 4530.640378074, 4510.298968946, 4479.400496614, 4433.605169893, 4367.694523391, 4275.729032416, 4151.075633815, 3986.010514103, 3770.536770139, 3490.245977189, 3123.800206228, 2642.327425121, 2016.860336457, 1251.15999803, 471.9330526949, 15.72497414425, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
varlist['TI_equil']=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 419.3295522232, 1399.265659004, 2318.863506893, 3010.402649893, 3500.700945034, 3850.597354456, 4103.298319454, 4283.677295642, 4407.319425897, 4486.677753339, 4533.381483685, 4558.298275196, 4570.466352163, 4575.803433794, 4576.953953376, 4574.166220348, 4565.04416692, 4541.399347934, 4484.198416614, 4363.27065221, 4146.049025335, 3803.274721578, 3297.424723466, 2570.704607261, 1583.83503038, 480.4338273235, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
varlist['BALLOON_grate']=np.array([-0.04038799836603, -0.04594611202707, -0.05159264035789, -0.05338423360233, -0.05538524916278, -0.06223263550467, -0.06593416618875, -0.06803765610926, -0.07049825141982, -0.07353553881503, -0.07627930828601, -0.07823856581361, -0.08028355417528, -0.0847518683282, -0.09088256686384, -0.096695078757, -0.1010977773673, -0.1042979533192, -0.1064567057502, -0.1079072692418, -0.1088968041487, -0.109626036973, -0.1103044820413, -0.1108589250405, -0.1114842619303, -0.1120915177094, -0.1124963882457, -0.1130936926979, -0.114902185858, -0.1185582143647, -0.1237070851507, -0.129403595285, -0.1349031265389, -0.1398443012642, -0.1441579758086, -0.1478751183043, -0.1509725891843, -0.1534357811374, -0.1553311534596, -0.1566305780363, -0.157322213026, -0.1572450696133, -0.1564453462638, -0.1546669915085, -0.152015521557, -0.1483808876154, -0.1441790729157, -0.1401886113428, -0.1381082648658, -0.1393894821852, -0.1439117986092, -0.1497296867413, -0.1544865746525, -0.157323965285, -0.1586332667534, -0.1594989383624, -0.1603046246896, -0.1614807744339, -0.1625760618772, -0.1639025855745, -0.1647907823608, -0.1655497099594, -0.1656995607924, -0.1656149533598, -0.1648199544707, -0.1635848840859, -0.1615789110373, -0.1588423550452, -0.1551415599088, -0.1502121596218, -0.1438606704745, -0.1362970859056, -0.1277817104033, -0.1192946351605, -0.1122512484435, -0.1079742988115, -0.1068304840944, -0.1079370212531, -0.1101260006045, -0.1120586151679, -0.1128430860159, -0.1115773140381, -0.107599709037, -0.1003344521832, -0.08995339321328, -0.08342165529579, -0.08473805102768, -0.08933490049766, -0.08905534025652, -0.07140087592814, -0.05655114596739, -0.06314969066222, -0.08752073467681, -0.1026093760905, -0.1115637043709, -0.1185383063478, -0.1242974403219, 0.0])
varlist['BOOTSTRAP_equil']=np.array([3741.616994374, 15025.80617273, 19428.83882039, 25237.13400923, 31039.81517359, 35667.86325406, 29972.24994714, 105262.6647243, 82478.40368427, 86887.58221936, 94463.97067028, 102872.2839846, 104824.0402004, 124808.539877, 133281.8878384, 130601.3284812, 208449.405322, 182241.7171977, 189790.5572342, 195325.7984364, 208492.0956995, 212649.9558268, 234490.3859978, 233328.1138398, 241301.6732594, 237736.0445714, 214870.1631477, 179275.4434528, 467.8708070572, 635561.8408471, 435152.6592574, 394574.7588066, 379888.7767066, 379510.5884794, 376680.8510696, 395123.4185516, 385816.5624314, 390165.7260593, 374378.6020315, 459659.3750827, 416377.6130999, 433431.4773128, 411500.6455102, 395831.5752467, 313217.9062821, 479938.9476415, 389431.9088218, 395432.7822432, 368713.6971767, 344698.998023, 408270.1324638, 333494.9037634, 251039.3828837, 284556.4823612, 171306.5370011, 130365.7954848, 39814.81495194, 87749.91776805, 302154.5929236, 712441.6132734, 1739215.436625, 7094384.578059, 4020579.822206, 2099490.400311, 1503628.586513, 1221716.568258, 1044189.972168, 934696.4368173, 864939.3426652, 723701.9637882, 712075.5929374, 595027.0088584, 681231.678732, 565657.1139803, 499845.9559692, 298795.5982667, 911023.4445974, 635725.7836139, 541297.9180987, 494254.7753182, 506836.9470338, 475673.5980304, 432666.1710154, 438943.4018987, 411202.3805122, 398234.1092202, 380641.2380251, 360641.9304723, 334209.7025736, 301458.4180866, 261857.7945433, 220771.9489361, 178851.3089605, 138655.4268259, 101540.8742553, 42946.30951772, 22.56759839734, 0.0])
varlist['NEO_equil']=np.array([2.559979695921e-07, 1.291925850969e-07, 1.56444233477e-07, 2.393404861728e-07, 3.45778820536e-07, 4.830332684432e-07, 6.549312025591e-07, 8.78936366412e-07, 1.113860755601e-06, 1.415732866551e-06, 1.728612404204e-06, 2.113944351495e-06, 2.572397689976e-06, 3.095179724915e-06, 3.615246396866e-06, 4.210078922015e-06, 4.890768993338e-06, 5.662054641078e-06, 6.657113572191e-06, 7.510249534397e-06, 8.507488068543e-06, 9.579235260336e-06, 1.074792306811e-05, 1.209537390499e-05, 1.332633994049e-05, 1.515165160415e-05, 1.67960193765e-05, 1.846885801059e-05, 2.077290484326e-05, 2.271790052823e-05, 2.524970974078e-05, 2.765630908478e-05, 3.075107121975e-05, 3.35925578673e-05, 3.68763422003e-05, 3.967599444416e-05, 4.419776510975e-05, 4.794474358461e-05, 5.255956025984e-05, 5.67718244444e-05, 6.150044764429e-05, 6.6621927914e-05, 7.233049824387e-05, 7.824584887813e-05, 8.445826910772e-05, 9.104636570827e-05, 9.780751517636e-05, 0.0001052838540742, 0.0001126136560128, 0.0001221159636624, 0.0001304256973182, 0.0001412951666099, 0.0001506801533267, 0.0001606720980109, 0.0001739951149556, 0.0001850059808813, 0.0001970664552271, 0.0002118436217842, 0.0002254716167018, 0.0002399312434909, 0.0002580940083635, 0.0002725163028699, 0.0002896509257039, 0.000306407412253, 0.0003247090403534, 0.0003433978408855, 0.0003629786734248, 0.0003873149865173, 0.0004086032598368, 0.0004331215090813, 0.0004572196347519, 0.0004820694781347, 0.0005091056984443, 0.0005395944216874, 0.0005692001397824, 0.0005979440596762, 0.0006341842969603, 0.0006757623219081, 0.0007160631678786, 0.0007537121955172, 0.0008034933002606, 0.0008511945194389, 0.0008992537638418, 0.0009592730767802, 0.001016381772563, 0.001079542967339, 0.001150939852542, 0.001224758617235, 0.00128835017095, 0.001371831338544, 0.001470435921633, 0.001557556612008, 0.001651103481077, 0.001761077542858, 0.001848772910163, 0.001969610689062, 0.002073455689253, 0.002208645719896])
varlist['TXPORT_equil']=np.array([0.09588515305729, 0.09653944559713, 0.09552027507684, 0.09543900059419, 0.09664499501116, 0.09807106048981, 0.09980965711076, 0.1011854860248, 0.1008589488946, 0.0963976919595, 0.09219282238409, 0.08853422431189, 0.08518680507368, 0.08578408063686, 0.09051342384874])

#print(data['HELICITY_FULL_equil'][25].tolist())
for temp in varlist:
    act = varlist[temp]
    cal = data[temp]
    cal = np.where(act==0,0,cal)
    div = np.where(act==0,1,act)
    perct = 100*sum(abs(act-cal)/div)
    print('  '+temp+': '+str(max(cal))+'   '+str(max(act))+'   '+str(int(perct))+'%')
    if perct > failtol:
        lfail = 1
print('=================')

if (lfail):
    print('  STATUS: FAIL!!!!!')
else:
    print('  STATUS: PASS')

sys.exit(0)



