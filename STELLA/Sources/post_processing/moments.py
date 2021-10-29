import numpy as np
from matplotlib import pyplot as plt
import stella_data as sd
from stella_plots import plot_1d
from matplotlib.backends.backend_pdf import PdfPages

for i in range(sd.nspec):
    plot_1d(sd.time,np.sqrt(sd.density[:,i,sd.iz0,0,0,1]**2+sd.density[:,i,sd.iz0,0,0,0]**2),xlab='time',ylab='density')
    plot_1d(sd.time,np.sqrt(sd.upar[:,i,sd.iz0,0,0,1]**2+sd.upar[:,i,sd.iz0,0,0,0]**2),xlab='time',ylab='upar')
    plot_1d(sd.time,np.sqrt(sd.temperature[:,i,sd.iz0,0,0,1]**2+sd.temperature[:,i,sd.iz0,0,0,0]**2)/np.sqrt(sd.temperature[0,i,sd.iz0,0,0,1]**2+sd.temperature[0,i,sd.iz0,0,0,0]**2),xlab='time',ylab='temperature')

file = sd.outdir+'moments.pdf'
pdf = PdfPages(file)
for i in plt.get_fignums():
    pdf.savefig(i)
pdf.close()
plt.close('all')
