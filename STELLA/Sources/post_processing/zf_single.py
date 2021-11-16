

from scipy.io import netcdf
from scipy.interpolate import UnivariateSpline
import numpy as np
import numpy.matlib

tave=1000
tfave=1090
fpref='center'

basedir = '/marconi_work/FUA34_MULTEI/stonge0_FUA34/rad_test/2nd_deriv/rho_scan2/r0.001ky/'
basedir = '/marconi_work/FUA34_MULTEI/stonge0_FUA34/rad_test/fg_drive/0.001d_alt/'
#basedir = '/marconi_work/FUA34_MULTEI/stonge0_FUA34/rad_test/fg_drive/0.001t_2/'
basedir = '/marconi_work/FUA34_MULTEI/stonge0_FUA34/rad_test/fg_drive/0.001d/'
#basedir = '/marconi_work/FUA34_MULTEI/stonge0_FUA34/rad_test/fg_drive/rhoscan_hr/0.002/'
#basedir = '/marconi_work/FUA35_OXGK/stonge0/rad_scan/0.001d_aht/'
basedir = '/marconi_work/FUA35_OXGK/stonge0/rad_scan/0.001d_per/'
center_file = basedir +fpref + '.out.nc'


center_nc = netcdf.netcdf_file(center_file,'r')


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

def phi_vs_t(infile,var,ny,nx):
# t ntube z kx ky ri

  avt, present = read_stella_float(infile,var)
  #print('c')
  arr = ny*nx*(avt[:,0,:,:,:,0] + 1j*avt[:,0,:,:,:,1])
  #print('d')
  #arr = np.fft.ifft(avt_kxky,axis=2)
  #print('e')

  return arr

naky  = center_nc.dimensions['ky']
nakxc = center_nc.dimensions['kx']

ky  = np.copy(center_nc.variables['ky'][:])
kxc  = np.copy(center_nc.variables['kx'][:])

Lxc = 2.*np.pi/kxc[1]
dx = 2*np.pi/kxc[1]/nakxc

t  = np.copy(center_nc.variables['t'][:])
nt = t.size

zed  = np.copy(center_nc.variables['zed'][:])
radgrid  = np.copy(center_nc.variables['rad_grid'][:])
kperp2 = np.copy(center_nc.variables['kperp2'][:])

print(str(np.amax(kperp2)) + ' ' + str(np.amax(ky)**2 + np.amax(kxc)**2))
nzed = zed.size
omp = ((nzed+1)/2) - 1
delzed = zed[1]-zed[0]
jacobc  = np.copy(center_nc.variables['jacob'][:])

dl_over_bc = np.squeeze(delzed*jacobc)
dl_over_bc[nzed-1] = 0.0
dl_over_bc = dl_over_bc/sum(dl_over_bc)

dobc = np.transpose(np.matlib.tile(dl_over_bc,(naky,nakxc,1)))

phic_kxky = phi_vs_t(center_nc,'phi_vs_t',naky,nakxc)
phiZF_kxky= np.sum(dobc[:,:,0]*phic_kxky[:,:,:,0],1)

#print(phic_kxky[nt-1,omp,:,0])
#print(phiZF_kxky[nt-1,:])

phizfc_k= np.sum((1j*kxc)**0*dobc[:,:,0]*phic_kxky[:,:,:,0],1)
vxc_k   = np.sum((1j*kxc)**1*dobc[:,:,0]*phic_kxky[:,:,:,0],1)
exbc_k  = np.sum((1j*kxc)**2*dobc[:,:,0]*phic_kxky[:,:,:,0],1)


phizfc = np.real(np.fft.ifft(phizfc_k,axis=1))/naky
vxc    = np.real(np.fft.ifft(vxc_k,axis=1))/naky
exbc   = np.real(np.fft.ifft(exbc_k,axis=1))/naky
exbsmooth = np.copy(exbc)
#vxc = np.fft.ifft(vxc_k,axis=2)

tind=nt-1
nti=nt
for i in range (0, nt):
  if(t[i]> tave):
    tind = i
    break
for i in range (0, nt):
  if(t[i]> tfave):
    nti = i
    break
p_ave   = np.mean(phizfc[tind:nti,:],0)
v_ave   = np.mean(vxc[tind:nti,:],0)
e_ave   = np.squeeze(np.mean(exbc[tind:nti,:],0))
spec_ave=np.mean(np.abs(exbc_k[tind:nti]),0)

window_width=7
w2 = (window_width+1)/2
cumsum_vec = numpy.cumsum(numpy.insert(e_ave, 0, 0))
ma_vec = (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width
spl = UnivariateSpline(radgrid[:,0],e_ave)
spl.set_smoothing_factor(0.5)
ma_vec = spl(radgrid[:,0])

ma_vec[0]      = 0.5*(e_ave[0] +e_ave[1]) 
ma_vec[nakxc-1] = 0.5*(e_ave[nakxc-2] +e_ave[nakxc-1]) 
for i in range (1,e_ave.size-1):
  ma_vec[i] = 0.25*e_ave[i-1]+0.5*e_ave[i]+0.25*e_ave[i+1]
  exbsmooth[:,i] =0.25*exbc[:,i-1]+0.5*exbc[:,i] + 0.25*exbc[:,i+1]
#for i in range (e_ave.size):
# ss=0
# sl=0
# count=0
# for j in range (max(0,i-1),min(e_ave.size,i+1)):
#     count=count+1
#     ss+=e_ave[j]
#     sl+=exbc[:,j]
# ma_vec[i] = ss/count  
# exbsmooth[:,i] = sl/count


cout = open(basedir + fpref + '.exb_t','w')
cout.write('[1] t     ')
cout.write('[2] x     ')
cout.write('[3] x simp')
cout.write('[4] rho   ')
cout.write('[5] phi   ')
cout.write('[6] v     ')
cout.write('[7] exb   ')
cout.write('[8] exb s ')
cout.write('[9] phik  ')
cout.write('\n')
for i in range (nt):
  for j in range (nakxc):
    cout.write('%e ' % t[i])
    cout.write('%e ' % radgrid[j,0])
    cout.write('%e ' % (dx*j))
    cout.write('%e ' % radgrid[j,1])
    cout.write('%e ' % phizfc[i,j])
    cout.write('%e ' % vxc[i,j])
    cout.write('%e ' % exbc[i,j])
    cout.write('%e ' % exbsmooth[i,j])
    cout.write('%e ' % abs(phizfc_k[i,j]))
    cout.write('\n')
  cout.write('\n')

cout.close()

cout = open(basedir + fpref + '.exb_ave','w')
cout.write('#Average from t=' + str(t[tind])+ ' to t=' + str(t[nti-1]) + '\n')
cout.write('[1] x     ')
cout.write('[2] phi   ')
cout.write('[3] v     ')
cout.write('[4] exb   ')
cout.write('[5] exb smooth')
cout.write('[6] spec ave')
cout.write('\n')
for j in range (w2,nakxc-w2):
  cout.write('%e ' % radgrid[j,0])
  cout.write('%e ' % p_ave[j])
  cout.write('%e ' % v_ave[j])
  cout.write('%e ' % e_ave[j])
  cout.write('%e ' % ma_vec[j])
  cout.write('%e ' % spec_ave[j])
  cout.write('\n')
cout.close()

cout = open(basedir + fpref + '.exb_en','w')
cout.write('[1] t     ')
cout.write('[2] phi2  ')
cout.write('[3] uzf2  ')
cout.write('[4] exb2  ')
cout.write('\n')
for i in range (nt):
    cout.write('%e ' % t[i])
    cout.write('%e ' % np.mean(phizfc[i,:]**2))
    cout.write('%e ' % np.mean(vxc[i,:]**2))
    cout.write('%e ' % np.mean(exbc[i,:]**2))
    cout.write('\n')
cout.close()

