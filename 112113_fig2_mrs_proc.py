#!/usr/bin/python
import matplotlib               
import numpy                    
from pylab import *
from matplotlib import pyplot as pl
from matplotlib import axis as ax
from matplotlib.ticker import MultipleLocator#, FormatStrFormatter
#from scipy import integrate
from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('plots/130805_Hopkins_y_with_ymax_inset.pdf')

x_nopeak = numpy.loadtxt('data/x_nopeak.txt', unpack=True, usecols = [0])
y_nopeak = numpy.loadtxt('data/a-sio2-eps2.txt', unpack=True, usecols = [1])

x_peak = numpy.loadtxt('data/x_peak.txt', unpack=True, usecols = [0])
y_peak = numpy.loadtxt('data/Y24204L.txt', unpack=True, usecols = [1])# L => eps2


eiz_nopeak = numpy.loadtxt('data/eiz_nopeak.txt', unpack=True, usecols = [0])
eiz_peak = numpy.loadtxt('data/eiz_peak.txt', unpack=True, usecols = [0])


r09,r10,r15,r20 = numpy.loadtxt('data/131128_inset_MRS_fig.txt', unpack=True, usecols = [0,1,2,3])

z_y = linspace(0,10,250)
sum_Li2_A_dep = numpy.empty(len(z_y))
yA_dep = numpy.empty(len(z_y))
norm_yA_dep = numpy.empty(len(z_y))
A_dep = [2.,4.,6.,8.,10.,12.]#linspace(5,36,6)
y_0 = numpy.empty(len(A_dep))
sum_dF =  0.0
sum_F =  0.0
sum_F_p =  0.0
diff_sum_F =  0.0
max_yA = numpy.empty(len(A_dep))
max_xA = numpy.empty(len(A_dep))

fig = pl.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8])
colors = ['b','g','r','c','m','y','k']
for r in range(len(A_dep)):
    labels=r'$A\, =\, %.2f$' % A_dep[r]
    for t in range(len(z_y)):
        sum_Li2_A_dep[t] = 0.0
        for u in arange(1,101):
            sum_Li2_A_dep[t] += (((A_dep[r]/(A_dep[r]+2.*z_y[t]**2 +2))**2)**u)/(u**2)
        yA_dep[t] = (sum_Li2_A_dep[t])*\
        ((1. - ((A_dep[r]/(1.+z_y[t]**2)) / (2.+(A_dep[r]/(1.+z_y[t]**2))))**2)/\
        ((A_dep[r]/(1.+z_y[t]**2)) / (2.+(A_dep[r]/(1.+z_y[t]**2)))))*\
        (1./(1.+z_y[t]**2))/(1. + A_dep[r]/(1+z_y[t]**2))
 
#        yA_dep[t] = (sum_Li2_A_dep[t]*\
#        (A_dep[r]/(A_dep[r]+2.*z_y[t]**2+2))**(-1)*\
#        (1-(A_dep[r]/(A_dep[r]+2.*z_y[t]**2+2))**2))\

#    max_xA[r]= z_y[yA_dep.argmax()]
    
    # MAIN PLOT
    # ----------------------------------------------------------
    ax.plot(z_y,yA_dep/yA_dep[0],  color = colors[r], label = labels)
    max_yA[r]= max(yA_dep)/yA_dep[0]	
    max_xA[r]= z_y[yA_dep.argmax()]
    ax.plot(max_xA, max_yA, linestyle = '', marker = 's', markerfacecolor='m')#,linestyle = lts[i],label=labels)#color = colors[j],

pl.xlabel(r'$\zeta/\omega_{0}$', size = 21)
pl.ylabel(r'$ \mathcal{\delta A}(i\zeta/\omega_{0} ) $', size = 21)
# INSET
# ----------------------------------------------------------
ax_inset = fig.add_axes([0.5,0.5,0.36,0.36])
#ax_inset = fig.add_axes([0.55,0.55,0.33,0.33])
# NB, if you change the above eqn for y, you must regenerate this input:
A_dep_r = linspace(0.0000000000001,10,300)#(0.000001,10,20)
#Ar = [A_dep_r,A_dep_r,A_dep_r,A_dep_r]
#r_dep = [r09,r10,r15,r20]

#for i in len(r_dep):
#ax_inset.plot(Ar[i],r_dep[i], color = 'k', linestyle = ':')
ax_inset.plot(A_dep_r,r09,'k:',A_dep_r,r10,'m-',A_dep_r,r15,'k--',A_dep_r,r20,'k--')
#ax_inset.plot(A_dep_r,r09,'k:',A_dep_r,r10,'m-',A_dep_r,r15,'k-.',A_dep_r,r20,'k--')
#ax_inset.plot(A_dep,r10, color = 'k', linestyle = '-.')
#ax_inset.plot(A_dep,r15, color = 'k', linestyle = '--')
#ax_inset.plot(A_dep,r20, color = 'k', linestyle = '-')
pl.xlabel(r'$\epsilon(0)-1$',size = 14)#\,=\,\frac{2}{\pi}C
pl.ylabel(r'$\rm{MAX}$ $[\zeta/\omega_{0}] $',size = 14)
pl.tick_params(labelsize = 'small')
#ax.legend(loc = 'best')
pl.savefig('plots/MRS_proc_fig2.pdf')

pp.savefig()
pl.show()
pp.close()

pl.close()        

