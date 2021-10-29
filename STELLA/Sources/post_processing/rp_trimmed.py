

from scipy.io import netcdf
import numpy as np
import numpy.matlib


tave = 450
br = 16

basedir = '/marconi_work/FUA34_MULTEI/stonge0_FUA34/rad_test/2nd_deriv/T1/'
basedir = '/marconi_work/FUA34_MULTEI/stonge0_FUA34/rad_test/2nd_deriv/rho_scan2/r0.001/'
basedir = '/marconi_work/FUA34_MULTEI/stonge0_FUA34/rad_test/fg_drive/rhoscan_hr/0.002/'
basedir = '/marconi_work/FUA34_MULTEI/stonge0_FUA34/rad_test/fg_drive/0.001d/'
basedir = '/marconi_work/FUA35_OXGK/stonge0/rad_scan/source_test/proj6/'
#basedir = '/marconi_work/FUA35_OXGK/stonge0/rad_scan/new_source_ky/'
basedir = '/marconi_work/FUA35_OXGK/stonge0/rad_scan/CBC_flux_alt/'
basedir = '/marconi_work/FUA35_OXGK/stonge0/rad_scan/new_source_hr_v3e_alt/'
basedir = '/marconi_work/FUA35_OXGK/stonge0/rad_scan/new_source_hr_v3e_pe3/'
basedir = '/marconi_work/FUA35_OXGK/stonge0/rad_scan/rho_scan/0.005/'
left_file   = basedir + 'left.out.nc'
center_file = basedir + 'center.out.nc'
right_file   = basedir + 'right.out.nc'


right_nc  = netcdf.netcdf_file(right_file,'r')
center_nc = netcdf.netcdf_file(center_file,'r')
left_nc   = netcdf.netcdf_file(left_file,'r')


def read_stella_float(infile, var):

  import numpy as np

  try:
    #print('a')
    #arr = np.copy(infile.variables[var][:])
    arr = infile.variables[var][:]
    #print('b')
    flag = True
  except KeyError:
    print('INFO: '+var+' not found in netcdf file')
    arr =np.arange(1,dtype=float)
    flag = FLAG

  return arr, flag

print('0')
naky  = center_nc.dimensions['ky']
nakxl =   left_nc.dimensions['kx']
nakxc = center_nc.dimensions['kx']
nakxr =  right_nc.dimensions['kx']
ky  = np.copy(center_nc.variables['ky'][:])
kxc = np.copy(center_nc.variables['kx'][:])

t  = np.copy(center_nc.variables['t'][:])
nt = t.size

Lxc = 2.*np.pi/kxc[1]
dxc = Lxc / nakxc

zed  = np.copy(center_nc.variables['zed'][:])
nzed = zed.size
omp = ((nzed+1)/2) - 1
delzed = zed[1]-zed[0]

radgrid  = np.copy(center_nc.variables['rad_grid'][:])

fac = 2*np.ones(naky)
fac[0] = 1

jacobl  = np.copy(  left_nc.variables['jacob'][:])
jacobc  = np.copy(center_nc.variables['jacob'][:])
jacobr  = np.copy( right_nc.variables['jacob'][:])

print('1')

dl_over_bl = np.squeeze(delzed*jacobl)
dl_over_bc = np.squeeze(delzed*jacobc)
dl_over_br = np.squeeze(delzed*jacobr)
dl_over_bl[nzed-1] = 0.0
dl_over_bc[nzed-1] = 0.0
dl_over_br[nzed-1] = 0.0
dl_over_bl = dl_over_bl/sum(dl_over_bl)
dl_over_bc = dl_over_bc/sum(dl_over_bc)
dl_over_br = dl_over_br/sum(dl_over_br)

dobl = np.transpose(np.matlib.tile(dl_over_bl,(naky,nakxl,1)))
dobc = np.transpose(np.matlib.tile(dl_over_bc,(naky,nakxc,1)))
dobr = np.transpose(np.matlib.tile(dl_over_br,(naky,nakxr,1)))

print('2')

# t spec x
pfluxl  = np.copy(  left_nc.variables['pflux_x'][:])
pfluxc  = np.copy(center_nc.variables['pflux_x'][:])
pfluxr  = np.copy( right_nc.variables['pflux_x'][:])

vfluxl  = np.copy(  left_nc.variables['vflux_x'][:])
vfluxc  = np.copy(center_nc.variables['vflux_x'][:])
vfluxr  = np.copy( right_nc.variables['vflux_x'][:])

qfluxl  = np.copy(  left_nc.variables['qflux_x'][:])
qfluxc  = np.copy(center_nc.variables['qflux_x'][:])
qfluxr  = np.copy( right_nc.variables['qflux_x'][:])

print('3')

densl  = np.copy(left_nc.variables['dens_x'][:])
uparl  = np.copy(left_nc.variables['upar_x'][:])
templ  = np.copy(left_nc.variables['temp_x'][:])

densc  = np.copy(center_nc.variables['dens_x'][:])
uparc  = np.copy(center_nc.variables['upar_x'][:])
tempc  = np.copy(center_nc.variables['temp_x'][:])

densr  = np.copy(right_nc.variables['dens_x'][:])
uparr  = np.copy(right_nc.variables['upar_x'][:])
tempr  = np.copy(right_nc.variables['temp_x'][:])

cout = open(basedir + 'left.fluxes_t','w')
cout.write('[1] t     ')
cout.write('[2] x     ')
cout.write('[3] flux_d')
cout.write('[4] flux_u')
cout.write('[5] flux_t')
cout.write('\n')
print('5')
for i in range (0, nt):
  for j in range (0, nakxl):
    cout.write('%e ' % t[i])
    cout.write('%e ' % (dxc*j))
    cout.write('%e ' % pfluxl[i,0,j])
    cout.write('%e ' % vfluxl[i,0,j])
    cout.write('%e ' % qfluxl[i,0,j])
    cout.write('\n')
  cout.write('\n')
cout.close()

cout = open(basedir + 'center.fluxes_t','w')
cout.write('[1] t     ')
cout.write('[2] x     ')
cout.write('[3] x simp')
cout.write('[4] r     ')
cout.write('[5] flux_d')
cout.write('[6] flux_u')
cout.write('[7] flux_t')
cout.write('\n')
print('6')
for i in range (0, nt):
  for j in range (0, nakxc):
    cout.write('%e ' % t[i])
    cout.write('%e ' % radgrid[j,0]) 
    cout.write('%e ' % (dxc*j))
    cout.write('%e ' % radgrid[j,2]) 
    cout.write('%e ' % pfluxc[i,0,j])
    cout.write('%e ' % vfluxc[i,0,j])
    cout.write('%e ' % qfluxc[i,0,j])
    cout.write('%e ' % densc[i,0,j])
    cout.write('%e ' % uparc[i,0,j])
    cout.write('%e ' % tempc[i,0,j])
    cout.write('\n')
  cout.write('\n')
cout.close()

for i in range (0, nt):
  cout = open(basedir + 'center.fluxes_' + str(i),'w')
  for j in range (0, nakxc):
    cout.write('%e ' % radgrid[j,2]) 
    cout.write('%e ' % tempc[i,0,j])
    cout.write('\n')
  cout.close()

cout = open(basedir + 'right.fluxes_t','w')
cout.write('[1] t     ')
cout.write('[2] x     ')
cout.write('[3] flux_d')
cout.write('[4] flux_u')
cout.write('[5] flux_t')
cout.write('\n')
print('7')
for i in range (0, nt):
  for j in range (0, nakxr):
    cout.write('%e ' % t[i])
    cout.write('%e ' % (dxc*j))
    cout.write('%e ' % pfluxr[i,0,j])
    cout.write('%e ' % vfluxr[i,0,j])
    cout.write('%e ' % qfluxr[i,0,j])
    cout.write('\n')
  cout.write('\n')
cout.close()

tind=nt-1
for i in range (0, nt):
  if(t[i]> tave):
    tind = i
    break
    
print(str(tind) + '  ' + str(nt))

print('8')

plave = np.mean(pfluxl[tind:nt,0,:],0)
vlave = np.mean(vfluxl[tind:nt,0,:],0)
qlave = np.mean(qfluxl[tind:nt,0,:],0)

pcave = np.mean(pfluxc[tind:nt,0,:],0)
vcave = np.mean(vfluxc[tind:nt,0,:],0)
qcave = np.mean(qfluxc[tind:nt,0,:],0)

prave = np.mean(pfluxr[tind:nt,0,:],0)
vrave = np.mean(vfluxr[tind:nt,0,:],0)
qrave = np.mean(qfluxr[tind:nt,0,:],0)

dlave = np.mean(densl[tind:nt,0,:],0)
ulave = np.mean(uparl[tind:nt,0,:],0)
tlave = np.mean(templ[tind:nt,0,:],0)

dcave = np.mean(densc[tind:nt,0,:],0)
ucave = np.mean(uparc[tind:nt,0,:],0)
tcave = np.mean(tempc[tind:nt,0,:],0)

drave = np.mean(densr[tind:nt,0,:],0)
urave = np.mean(uparr[tind:nt,0,:],0)
trave = np.mean(tempr[tind:nt,0,:],0)

print('9')

dens_zero = np.mean(densc[:,0,:],1)
upar_zero = np.mean(uparc[:,0,:],1)
temp_zero = np.mean(tempc[:,0,:],1)

dens_br_zero = np.mean(densc[:,0,br:(nakxc-br-1)],1)
upar_br_zero = np.mean(uparc[:,0,br:(nakxc-br-1)],1)
temp_br_zero = np.mean(tempc[:,0,br:(nakxc-br-1)],1)

deri1 = np.arange(1,nakxc-2*br,1)
deri1 = deri1 - np.mean(deri1)
deri1 = deri1/np.sqrt(np.sum(deri1**2))

dens_br_d1 = np.sum(deri1*densc[:,0,br:(nakxc-br-1)],1)
upar_br_d1 = np.sum(deri1*uparc[:,0,br:(nakxc-br-1)],1)
temp_br_d1 = np.sum(deri1*tempc[:,0,br:(nakxc-br-1)],1)

cout = open(basedir + 'center.prof_ave','w')
cout.write('#Average from t=' + str(t[tind])+ ' to t=' + str(t[nt-1]) + '\n')
cout.write('#')
cout.write('[1] x       ')
cout.write('[2] x simp  ')
cout.write('[3] r       ')
cout.write('[4] pflux   ')
cout.write('[5] vflux   ')
cout.write('[6] qflux   ')
cout.write('[7] dens    ')
cout.write('[8] upar    ')
cout.write('[9] temp    ')
cout.write('\n')

for i in range (0, nakxc):
  cout.write('%e ' % radgrid[i,0])
  cout.write('%e ' % (dxc*i))
  cout.write('%e ' % radgrid[i,2])
  cout.write('%e ' % pcave[i])
  cout.write('%e ' % vcave[i])
  cout.write('%e ' % qcave[i])
  cout.write('%e ' % dcave[i])
  cout.write('%e ' % ucave[i])
  cout.write('%e ' % tcave[i])
  cout.write('\n')

cout.close()
cout = open(basedir + 'left.prof_ave','w')
cout.write('#Average from t=' + str(t[tind])+ ' to t=' + str(t[nt-1]) + '\n')
cout.write('#')
cout.write('[1] x       ')
cout.write('[2] pflux   ')
cout.write('[3] vflux   ')
cout.write('[4] qflux   ')
cout.write('\n')

for i in range (0, nakxl):
  cout.write('%e ' % (dxc*i-0.5*Lxc))
  cout.write('%e ' % plave[i])
  cout.write('%e ' % vlave[i])
  cout.write('%e ' % qlave[i])
  cout.write('%e ' % dlave[i])
  cout.write('%e ' % ulave[i])
  cout.write('%e ' % tlave[i])
  cout.write('\n')

cout.close()
cout.close()
cout = open(basedir + 'right.prof_ave','w')
cout.write('#Average from t=' + str(t[tind])+ ' to t=' + str(t[nt-1]) + '\n')
cout.write('#')
cout.write('[1] x       ')
cout.write('[2] pflux   ')
cout.write('[3] vflux   ')
cout.write('[4] qflux   ')
cout.write('\n')

for i in range (0, nakxr):
  cout.write('%e ' % (dxc*i-0.5*Lxc))
  cout.write('%e ' % prave[i])
  cout.write('%e ' % vrave[i])
  cout.write('%e ' % qrave[i])
  cout.write('%e ' % drave[i])
  cout.write('%e ' % urave[i])
  cout.write('%e ' % trave[i])
  cout.write('\n')

cout.close()

cout = open(basedir + 'center.zero_mode','w')
cout.write('[1] t    ')
cout.write('[2] dens ')
cout.write('[3] upar ')
cout.write('[4] temp ')
cout.write('\n')
print(basedir + 'center.zero_mode')
for i in range (0, nt):
  cout.write('%e ' % t[i])
  cout.write('%e ' % dens_zero[i]) 
  cout.write('%e ' % upar_zero[i])
  cout.write('%e ' % temp_zero[i]) 
  cout.write('%e ' % dens_br_zero[i]) 
  cout.write('%e ' % upar_br_zero[i])
  cout.write('%e ' % temp_br_zero[i]) 
  cout.write('%e ' % dens_br_d1[i]) 
  cout.write('%e ' % upar_br_d1[i])
  cout.write('%e ' % temp_br_d1[i]) 
  cout.write('\n')
cout.close()
