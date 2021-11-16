

from scipy.io import netcdf
import numpy as np
import numpy.matlib


basedir = '/marconi_scratch/userexternal/dstonge0/stella/rad_test/krook/krook1c_alt/'
basedir = '/marconi_work/FUA34_MULTEI/stonge0_FUA34/rad_test/2nd_deriv/rho_scan2/r0.001ky/'
basedir = '/marconi_work/FUA34_MULTEI/stonge0_FUA34/rad_test/fg_drive/0.001d/'
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

def phi_vs_t_to_x(infile,var,ny,nx):
# t ntube z kx ky ri

  avt, present = read_stella_float(infile,var)
  #print('c')
  avt_kxky = ny*nx*(avt[:,0,:,:,:,0] + 1j*avt[:,0,:,:,:,1])
  #print('d')
  arr = np.fft.ifft(avt_kxky,axis=2)
  #print('e')

  return arr

def mom_vs_t_to_x(infile,var,ny,nx):
  #in:  t nspec ntube z kx ky ri
  #out: t z kx ky

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
kxc  = np.copy(center_nc.variables['kx'][:])

dx = 2*np.pi/kxc[1]/nakxc

t  = np.copy(center_nc.variables['t'][:])
nt = t.size

zed  = np.copy(center_nc.variables['zed'][:])
nzed = zed.size
omp = ((nzed+1)/2) - 1
delzed = zed[1]-zed[0]

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

phil_xky = phi_vs_t_to_x(left_nc  ,'phi_vs_t',naky,nakxl)
phic_xky = phi_vs_t_to_x(center_nc,'phi_vs_t',naky,nakxc)
phir_xky = phi_vs_t_to_x(right_nc ,'phi_vs_t',naky,nakxr)

print('3')

densl_xky = mom_vs_t_to_x(left_nc  ,'density',naky,nakxl)
densc_xky = mom_vs_t_to_x(center_nc,'density',naky,nakxc)
densr_xky = mom_vs_t_to_x(right_nc ,'density',naky,nakxr)


uparl_xky = mom_vs_t_to_x(left_nc  ,'upar',naky,nakxl)
uparc_xky = mom_vs_t_to_x(center_nc,'upar',naky,nakxc)
uparr_xky = mom_vs_t_to_x(right_nc ,'upar',naky,nakxr)

templ_xky = mom_vs_t_to_x(left_nc  ,'temperature',naky,nakxl)
tempc_xky = mom_vs_t_to_x(center_nc,'temperature',naky,nakxc)
tempr_xky = mom_vs_t_to_x(right_nc ,'temperature',naky,nakxr)

vxl = 1j*ky*phil_xky
vxc = 1j*ky*phic_xky
vxr = 1j*ky*phir_xky


print('4')
phic_zf = np.real(np.sum(dobc[:,:,0]*phic_xky[:,:,:,0],1))
dens_zf = np.real(np.sum(dobc[:,:,0]*densc_xky[:,:,:,0],1))
upar_zf = np.real(np.sum(dobc[:,:,0]*uparc_xky[:,:,:,0],1))
temp_zf = np.real(np.sum(dobc[:,:,0]*tempc_xky[:,:,:,0],1))

phic2 = np.real(np.mean(np.abs(phic_xky[:,omp,:,:])**2,2))
dens2 = np.real(np.mean(np.abs(densc_xky[:,omp,:,:])**2,2))
upar2 = np.real(np.mean(np.abs(uparc_xky[:,omp,:,:])**2,2))
temp2 = np.real(np.mean(np.abs(tempc_xky[:,omp,:,:])**2,2))

print('5')
fluxl_d = 0.5*np.mean(np.sum(fac*np.real(vxl*np.conj(densl_xky)*dobl),1),2)/naky
fluxl_u = 0.5*np.mean(np.sum(fac*np.real(vxl*np.conj(uparl_xky)*dobl),1),2)/naky
fluxl_T = 0.5*np.mean(np.sum(fac*np.real(vxl*np.conj(templ_xky)*dobl),1),2)/naky
cout = open(basedir + 'left.fluxes_t','w')
cout.write('[1] t     ')
cout.write('[2] x     ')
cout.write('[3] flux_d')
cout.write('[4] flux_u')
cout.write('[5] flux_t')
cout.write('\n')
print('6')
for i in range (0, nt):
  for j in range (0, nakxl):
    cout.write('%f ' % t[i])
    cout.write('%f ' % (dx*j))
    cout.write('%f ' % fluxl_d[i,j])
    cout.write('%f ' % fluxl_u[i,j])
    cout.write('%f ' % fluxl_T[i,j])
    cout.write('\n')
  cout.write('\n')
cout.close()

fluxc_d = 0.5*np.mean(np.sum(fac*np.real(vxc*np.conj(densc_xky)*dobc),1),2)/naky
fluxc_u = 0.5*np.mean(np.sum(fac*np.real(vxc*np.conj(uparc_xky)*dobc),1),2)/naky
fluxc_T = 0.5*np.mean(np.sum(fac*np.real(vxc*np.conj(tempc_xky)*dobc),1),2)/naky
cout = open(basedir + 'center.fluxes_t','w')
cout.write('[1] t     ')
cout.write('[2] x     ')
cout.write('[3] flux_d')
cout.write('[4] flux_u')
cout.write('[5] flux_t')
cout.write('\n')
print('7')
for i in range (0, nt):
  for j in range (0, nakxc):
    cout.write('%f ' % t[i])
    cout.write('%f ' % (dx*j))
    cout.write('%f ' % fluxc_d[i,j])
    cout.write('%f ' % fluxc_u[i,j])
    cout.write('%f ' % fluxc_T[i,j])
    cout.write('\n')
  cout.write('\n')
cout.close()

fluxr_d = 0.5*np.mean(np.sum(fac*np.real(vxr*np.conj(densr_xky)*dobr),1),2)/naky
fluxr_u = 0.5*np.mean(np.sum(fac*np.real(vxr*np.conj(uparr_xky)*dobr),1),2)/naky
fluxr_T = 0.5*np.mean(np.sum(fac*np.real(vxr*np.conj(tempr_xky)*dobr),1),2)/naky
cout = open(basedir + 'right.fluxes_t','w')
cout.write('[1] t     ')
cout.write('[2] x     ')
cout.write('[3] flux_d')
cout.write('[4] flux_u')
cout.write('[5] flux_t')
cout.write('\n')
print('8')
for i in range (0, nt):
  for j in range (0, nakxr):
    cout.write('%f ' % t[i])
    cout.write('%f ' % (j*dx))
    cout.write('%f ' % fluxr_d[i,j])
    cout.write('%f ' % fluxr_u[i,j])
    cout.write('%f ' % fluxr_T[i,j])
    cout.write('\n')
  cout.write('\n')
cout.close()


tave = int(0.7*float(nt))
print('Average time ' + str(tave) + ' to ' + str(nt))
temp_ave = np.mean(temp_zf[tave:nt,:],0)
fdc_ave  = np.mean(fluxc_d[tave:nt,:],0)
fuc_ave  = np.mean(fluxc_u[tave:nt,:],0)
fTc_ave  = np.mean(fluxc_T[tave:nt,:],0)


cout = open(basedir + 'center.stuff','w')
cout.write('[1] i       ')
cout.write('[2] phi     ')
cout.write('[3] dens    ')
cout.write('[4] upar    ')
cout.write('[5] temp    ')
cout.write('[6] temp_ave')
cout.write('[7] flux_d  ')
cout.write('[8] flux_u  ')
cout.write('[9] flux_T  ')
cout.write('[10] fd_ave ')
cout.write('[11] fu_ave ')
cout.write('[12] fT_ave ')
cout.write('\n')
print('9')
for i in range (0, nakxc):
  cout.write('%f ' % (i*dx))
  cout.write('%f ' % phic_zf[nt-1,i])
  cout.write('%f ' % dens_zf[nt-1,i])
  cout.write('%f ' % upar_zf[nt-1,i])
  cout.write('%f ' % temp_zf[nt-1,i])
  cout.write('%f ' % temp_ave[i])
  cout.write('%f ' % fluxc_d[nt-1,i])
  cout.write('%f ' % fluxc_u[nt-1,i])
  cout.write('%f ' % fluxc_T[nt-1,i])
  cout.write('%f ' % fdc_ave[i])
  cout.write('%f ' % fuc_ave[i])
  cout.write('%f ' % fTc_ave[i])
  cout.write('\n')
#  print(("%d " % i), end='')
#  print( "%f " % dens_zf[nt-1,i], end='', file=cout )
#  print( "%f " % upar_zf[nt-1,i], end='', file=cout )
#  print( "%f " % temp_zf[nt-1,i], end='', file=cout )
#  print( "%f " % temp_ave[i]    , end='', file=cout )
#  print( "" ,file=cout )
cout.close()

for j in range (0, nt):
  cout = open(basedir + 'prof_' + str(j),'w')
  for i in range (0, nakxc):
    cout.write('%e ' % (i*dx))
    cout.write('%f ' % phic_zf[j,i])
    cout.write('%f ' % dens_zf[j,i])
    cout.write('%f ' % upar_zf[j,i])
    cout.write('%f ' % temp_zf[j,i])
    cout.write('%f ' % phic2[j,i])
    cout.write('%f ' % dens2[j,i])
    cout.write('%f ' % upar2[j,i])
    cout.write('%f ' % temp2[j,i])
    cout.write('\n')
  cout.close()

temp0 = np.mean(temp_zf,1)

cout = open(basedir + 'temp.prof','w')
for i in range (0, nt - 1):
  cout.write('%e ' % t[i])
  cout.write('%e ' % temp0[i])
  cout.write('\n')
cout.close()
exit()
