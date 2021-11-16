import numpy as np
import matplotlib.pyplot as plt
from stella_data import zed, bmag, gradpar, outdir
from stella_plots import plot_1d
from matplotlib.backends.backend_pdf import PdfPages

plot_1d(zed,bmag,'z',ylab='bmag')
plot_1d(zed,gradpar,'z',ylab='gradpar')

file = outdir+'stella_geo.pdf'
pdf = PdfPages(file)
for i in plt.get_fignums():
    pdf.savefig(i)
pdf.close()
plt.close('all')
