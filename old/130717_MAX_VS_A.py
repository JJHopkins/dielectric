#!/usr/bin/python
import matplotlib               
import numpy                    
from pylab import *
from matplotlib import pyplot as pl
from matplotlib import axis as ax
from matplotlib.ticker import MultipleLocator#, FormatStrFormatter
from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('plots/130729_Hopkins_ymax_vs_A.pdf')

z_y = linspace(0,3,500)# z_y is z_n/w_0 ~= e16/e16
A_dep = linspace(0.001,10,400)
sum_Li2_A_dep = numpy.empty(len(z_y))
yA_dep = numpy.empty(len(z_y))

max_yA = numpy.empty(len(A_dep))
max_xA = numpy.empty(len(A_dep))

for r in range(len(A_dep)):
    for t in range(len(z_y)):
        sum_Li2_A_dep[t] = 0.0
        for u in arange(1,101):
            sum_Li2_A_dep[t] += (((A_dep[r]/(A_dep[r]+2.*z_y[t]**2 +2))**2)**u)/(u**2)
        
        yA_dep[t] = (sum_Li2_A_dep[t]*\
        (A_dep[r]/(A_dep[r]+2.*z_y[t]**2+2))**(-1)*\
        (1-(A_dep[r]/(A_dep[r]+2.*z_y[t]**2+2))**2))\

    max_xA[r]= z_y[yA_dep.argmax()]
data = zip(A_dep, max_xA)
savetxt('output/130729_Hopkins_ymax_output.txt', data)

pl.figure()
pl.plot(A_dep, max_xA)  
#pl.plot(A_dep,max_yA/yA_dep[0], color = colors[r])
pl.xlabel(r'$A\,=\,\frac{2}{\pi}C$', size = 'x-large')
pl.ylabel(r'$ \frac{\xi}{\omega_{0}}(y_{max}) $', size = 'x-large')
pl.axis([0, 10.2, 0.0, 2.0])
pp.savefig()
pl.show()
pl.close()

pp.close()
