## Description: This file reads in a .fluxes file and outputs important statistical
##              information, such as the average, standard deviation, and standard
##              error of the mean appropriate for a continuous time signal
##
##              This file assumes the data is equally spaced in time.

from scipy.io import netcdf
from scipy.interpolate import UnivariateSpline
import numpy as np
import numpy.matlib

def get_stats(data,cind_in):
    n_in = np.size(data)
    ave  = np.mean(data)
    std  = np.std(data)
    dfft  = np.fft.fft(data - ave,n=2*n_in)
    auto = np.real(np.fft.ifft(dfft*np.conj(dfft)))[0:n_in]
    auto = auto/auto[0]
    fac = 1 - np.linspace(0,n_in,num=n_in,endpoint=False)/n_in
    fac=auto*fac
    corr = 1 + 2*np.sum(fac[1:cind_in])
    sem = std/np.sqrt(n_in/corr)
    return ave, std, corr, sem, auto

#This is the time where we start averaging from
tave=200
tcorr=100

#open the file for reading
fpref='master'
basedir = '/Users/denisst-onge/stella/build/CBC_test/'
master_file = basedir +fpref + '.fluxes'
master_file = 'test.fluxes'
master_fluxes = np.transpose(np.loadtxt(master_file,skiprows=1))

time  = master_fluxes[0,:]
pflux = master_fluxes[1,:]
uflux = master_fluxes[2,:]
qflux = master_fluxes[3,:]

dt = time[1]-time[0]

#find the index where we want to start averaging our data
n = np.size(time)
tind=0
while (time[tind] < tave) and (tind < (n-1)):
    tind = tind + 1


psub=pflux[tind:]
usub=uflux[tind:] 
qsub=qflux[tind:]

subsize=np.size(psub)
cind=0
while (time[cind] < tcorr) and (cind < (subsize-1)):
    cind = cind + 1

pave, pstd, pcorr, psem, pauto = get_stats(psub,cind)
uave, ustd, ucorr, usem, uauto = get_stats(usub,cind)
qave, qstd, qcorr, qsem, qauto = get_stats(qsub,cind)

print('#pflux ave:' + str(pave) + ' stddev:' + str(pstd) + ' tcorr:' + str(pcorr*dt) +  ' SEM:' + str(psem))
print('#uflux ave:' + str(uave) + ' stddev:' + str(ustd) + ' tcorr:' + str(ucorr*dt) +  ' SEM:' + str(usem))
print('#qflux ave:' + str(qave) + ' stddev:' + str(qstd) + ' tcorr:' + str(qcorr*dt) +  ' SEM:' + str(qsem))
for i in range(subsize):
    print(str(time[i]) + ' ' + str(pauto[i]) + ' ' + str(uauto[i]) +  '  ' + str(qauto[i]))


