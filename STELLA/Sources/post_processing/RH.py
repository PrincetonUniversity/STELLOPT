

from scipy.io import netcdf
import numpy as np
import numpy.matlib

tave=20
basedir = '/marconi_work/FUA34_MULTEI/stonge0_FUA34/rad_test/RH/new_method/RH1_alt2/'
right_file  = basedir + 'center.out.nc'
center_file = basedir + 'center.out.nc'
left_file   = basedir + 'center.out.nc'


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

def phi_vs_t(infile,var,ny,nx):
# t ntube z kx ky ri

  avt, present = read_stella_float(infile,var)
  #print('c')
  arr = ny*nx*(avt[:,0,:,:,:,0] + 1j*avt[:,0,:,:,:,1])
  #print('d')
  arr = np.fft.ifft(arr,axis=2)
  #print('e')

  return arr

def mom_vs_t(infile,var,ny,nx):
# t nspec ntube z kx ky ri

  avt, present = read_stella_float(infile,var)
  avt_kxky = ny*nx*(avt[:,0,0,:,:,:,0] + 1j*avt[:,0,0,:,:,:,1])
  arr = np.fft.ifft(avt_kxky,axis=2)

  return arr

print('0')
naky  = center_nc.dimensions['ky']
nakxl =   left_nc.dimensions['kx']
nakxc = center_nc.dimensions['kx']
nakxr =  right_nc.dimensions['kx']
ky  = np.copy(center_nc.variables['ky'][:])

kxl  = np.copy(  left_nc.variables['kx'][:])
kxc  = np.copy(center_nc.variables['kx'][:])
kxr  = np.copy( right_nc.variables['kx'][:])

rad_gridc  = np.copy(center_nc.variables['rad_grid'][:])

Lxl = 2*np.pi/kxl[1]
Lxc = 2*np.pi/kxc[1]
Lxr = 2*np.pi/kxr[1]

print(Lxl, Lxc, Lxr)

dxl = Lxl/nakxl
dxc = Lxc/nakxc
dxr = Lxr/nakxr

print(dxl)

t  = np.copy(center_nc.variables['t'][:])
nt = t.size

zed  = np.copy(center_nc.variables['zed'][:])
nzed = zed.size
omp = ((nzed+1)/2) - 1
delzed = zed[1]-zed[0]

jacobl  = np.copy(  left_nc.variables['jacob'][:])
jacobc  = np.copy(center_nc.variables['jacob'][:])
jacobr  = np.copy( right_nc.variables['jacob'][:])

print(ky)

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

phil_xky = phi_vs_t(left_nc  ,'phi_vs_t',naky,nakxl)
phic_xky = phi_vs_t(center_nc,'phi_vs_t',naky,nakxc)
phir_xky = phi_vs_t(right_nc ,'phi_vs_t',naky,nakxr)
#phil_xky = mom_vs_t(left_nc  ,'density',naky,nakxl)
#phic_xky = mom_vs_t(center_nc,'density',naky,nakxc)
#phir_xky = mom_vs_t(right_nc ,'density',naky,nakxr)

phil_zf = np.real(np.sum(dobl[:,:,0]*phil_xky[:,:,:,0],1))
phic_zf = np.real(np.sum(dobc[:,:,0]*phic_xky[:,:,:,0],1))
phir_zf = np.real(np.sum(dobr[:,:,0]*phir_xky[:,:,:,0],1))

phil_ave = np.mean(phil_zf[(nt-1-tave):(nt-1),:],axis=0) 
phic_ave = np.mean(phic_zf[(nt-1-tave):(nt-1),:],axis=0) 
phir_ave = np.mean(phir_zf[(nt-1-tave):(nt-1),:],axis=0) 
print(phil_ave)

print(phic_zf[nt-1,:])

window_width=5
cumsum_vec = numpy.cumsum(numpy.insert(np.abs(phic_zf[0,:]), 0, 0)) 
ma_vec = (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width
print(ma_vec)
cumsum_vec = numpy.cumsum(numpy.insert(np.abs(phic_ave[:]), 0, 0)) 
ma_vec2 = (cumsum_vec[window_width:] - cumsum_vec[:-window_width]) / window_width
print(ma_vec2)



cout = open(basedir + 'left.ZF','w')
for i in range (0, nakxl):
  cout.write('%f ' % (dxl*i))
  cout.write('%f ' % phil_zf[0,i])
  cout.write('%f ' % phil_zf[nt-1,i])
  cout.write('\n')
cout.close()
cout = open(basedir + 'left.ZFt','w')
for i in range (0,nt):
  cout.write('%e ' % t[i])
  for j in range (0, nakxl):
    cout.write('%f ' % (phil_zf[i,j]/phil_zf[0,j]))
  cout.write('\n')
cout.close()

cout = open(basedir + 'center.ZF','w')
#for i in range (0, nakxc-window_width):
for i in range (0, nakxc):
  cout.write('%f ' % rad_gridc[i,0])
  cout.write('%f ' % rad_gridc[i,1])
  cout.write('%f ' % rad_gridc[i,2])
  cout.write('%f ' % phic_zf[0,i])
  cout.write('%f ' % phic_zf[nt-1,i])
  cout.write('%f ' % phic_ave[i])
#  cout.write('%f ' % ma_vec[i])
#  cout.write('%f ' % ma_vec2[i])
  cout.write('\n')
cout.close()
cout = open(basedir + 'center.ZFt','w')
for i in range (0,nt):
  cout.write('%e ' % t[i])
  for j in range (0, nakxc):
    cout.write('%f ' % (phic_zf[i,j]/phic_zf[0,j]))
  cout.write('\n')
cout.close()


cout = open(basedir + 'right.ZF','w')
for i in range (0, nakxr):
  cout.write('%f ' % (dxr*i))
  cout.write('%f ' % phir_zf[1,i])
  cout.write('%f ' % phir_zf[nt-1,i])
  cout.write('\n')
cout.close()
cout = open(basedir + 'right.ZFt','w')
for i in range (0,nt):
  cout.write('%e ' % t[i])
  for j in range (0, nakxr):
    cout.write('%f ' % (phir_zf[i,j]/phir_zf[0,j]))
  cout.write('\n')
cout.close()
