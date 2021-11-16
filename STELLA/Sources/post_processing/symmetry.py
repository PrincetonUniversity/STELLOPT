import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from stella_time import timeavg
from stella_plots import plot_2d
import stella_data as sd

def plot_gavg_vs_zvpa():

    cmap = 'YlGnBu'
    xlab = '$z$'
    ylab = '$v_{\parallel}$'
    title = 'avg $|g|^2$'

    gzvs_avg = np.arange(sd.zed.size*sd.vpa.size,dtype=float).reshape(sd.vpa.size,sd.zed.size)

    for i in range(sd.gzvs.shape[1]):
        for j in range(sd.gzvs.shape[2]):
            for k in range(sd.gzvs.shape[3]):
                gzvs_avg[j,k] = timeavg(sd.gzvs[:,i,j,k])*np.exp(2.*sd.vpa[j]**2)
        g_max = np.absolute(gzvs_avg).max()
        g_min = 0.0

        fig = plot_2d(gzvs_avg,sd.zed,sd.vpa,g_min,g_max,xlab,ylab,title+' (is= '+str(i+1)+')',cmap)

    return fig

plot_gavg_vs_zvpa()
file = sd.outdir+'symmetry.pdf'
pdf = PdfPages(file)
for i in plt.get_fignums():
    pdf.savefig(i)
pdf.close()
plt.close('all')
