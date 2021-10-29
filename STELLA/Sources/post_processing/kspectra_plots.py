import numpy as np
from matplotlib import pyplot as plt
from stella_data import time, ntime, outdir, nakx, naky, kx, ky
from matplotlib.backends.backend_pdf import PdfPages
from stella_data import phi2_vs_kxky as phi2
from stella_time import timeavg

phi2avg = np.arange(nakx*naky,dtype=float).reshape(nakx,naky)
fout=open(outdir+'kxky_spectra.txt','w+')
for ikx in range(phi2.shape[1]):
    for iky in range(phi2.shape[2]):
        phi2avg[ikx,iky] = timeavg(phi2[:,ikx,iky])
        fout.write(str(kx[ikx])+' '+str(ky[iky])+' '+str(phi2avg[ikx,iky])+'\n')
    fout.write('\n')
fout.close()

phi2ky = np.arange(naky,dtype=float)
f2out=open(outdir+'ky_spectra.txt','w+')
f2out.write('1) ky, 2) |phi|^2 \n')
for iky in range(naky):
    phi2ky[iky] = np.sum(phi2avg[:,iky])
    f2out.write(str(ky[iky])+' '+str(phi2ky[iky])+'\n')
ky0 = np.sum(phi2ky[1:]*ky[1:])/np.sum(phi2ky[1:])
ky1 = np.sum(phi2ky[1:])/np.sum(phi2ky[1:]/ky[1:])
f2out.write('#ky0: '+str(ky0)+' ky1: '+str(ky1)+'\n')
f2out.close()

fig=plt.figure(figsize=(12,8))

# distinguish between zonal/nonzonal modes
for ikx in range(phi2.shape[1]-1):
    plt.semilogy(time,phi2[:,ikx+1,0],'--')
    
for iky in range(phi2.shape[2]-1):
    for ikx in range(phi2.shape[1]):
        plt.semilogy(time,phi2[:,ikx,iky+1])
        
        plt.xlabel('$t (v_{t}/a)$')
        plt.xlim([time[0],time[ntime-1]])
        plt.title('$\Phi^2(k_x,k_y)$')

file = outdir+'stella_kspectra.pdf'
pdf = PdfPages(file)
for i in plt.get_fignums():
    pdf.savefig(i)
pdf.close()
plt.close('all')


