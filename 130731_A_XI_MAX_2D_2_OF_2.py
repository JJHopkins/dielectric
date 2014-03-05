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


x_ymax, y_ymax = numpy.loadtxt('output/130729_Hopkins_ymax_output.txt', unpack=True, usecols = [0,1])

#coeff = (2.41*1e14) # in rad/s
#n = arange(0,250)
#z = n * coeff
#A = 4.0 #[1,4,6]
#B = 7.0 #[1,4,6]
#C = 10.0 #[1,4,6]
#sum_Li2 = numpy.empty(len(z))
#dF = numpy.empty(len(z))
#sum_Li2_yA = numpy.empty(len(z_y))
#sum_Li2_yB = numpy.empty(len(z_y))
#sum_Li2_yC = numpy.empty(len(z_y))

#yA = numpy.empty(len(z_y))
#yB = numpy.empty(len(z_y))
#yC = numpy.empty(len(z_y))
#A_dep = linspace(1.,11.,300)

#sum_Li3 = numpy.empty(len(z))
#sum_Li3_p = numpy.empty(len(z))
#F_3 = numpy.empty(len(z))
#F_3_p = numpy.empty(len(z))

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

#for p in range(len(z_y)):
#    sum_Li2_yA[p] = 0.0
#    sum_Li2_yB[p] = 0.0
#    sum_Li2_yC[p] = 0.0
#    for q in arange(1,101):
#        #sum_Li2_yA[p] += (((A/(1+z_[p]))**q)/(q**2)# +  (((a/(a+2.*z_y[p]**2 +2)))**q)/(q**2))
#        sum_Li2_yA[p] += (((A/(A+2.*z_y[p]**2 +2))**2)**q)/(q**2)
#        sum_Li2_yB[p] += (((B/(B+2.*z_y[p]**2 +2))**2)**q)/(q**2)
#        sum_Li2_yC[p] += (((C/(C+2.*z_y[p]**2 +2))**2)**q)/(q**2)
#    yA[p] = sum_Li2_yA[p]*((A/(A+2.*z_y[p]**2 +2))**(-1)*(1 - (A/(A+2.*z_y[p]**2 +2))**2))   
#    yB[p] = sum_Li2_yB[p]*((B/(B+2.*z_y[p]**2 +2))**(-1)*(1 - (B/(B+2.*z_y[p]**2 +2))**2))  
#    yC[p] = sum_Li2_yC[p]*((C/(C+2.*z_y[p]**2 +2))**(-1)*(1 - (C/(C+2.*z_y[p]**2 +2))**2)) 


fig = pl.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8])
colors = ['b','c','g','y','r','k']
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
    max_yA[r]= max(yA_dep)/yA_dep[0]	
    max_xA[r]= z_y[yA_dep.argmax()]
    ax.plot(max_xA, max_yA, linestyle = '', marker = 's', markerfacecolor='m')#,linestyle = lts[i],label=labels)#color = colors[j],
    
    # MAIN PLOT
    # ----------------------------------------------------------
    ax.plot(z_y,yA_dep/yA_dep[0], linestyle = '-',  color = colors[r], label = labels)

pl.xlabel(r'$\zeta/\omega_{0}$', size = 21)
pl.ylabel(r'$ \mathcal{\delta A}(i\zeta/\omega_{0} ) $', size = 21)
# INSET
# ----------------------------------------------------------
ax_inset = fig.add_axes([0.5,0.5,0.36,0.36])
#ax_inset = fig.add_axes([0.55,0.55,0.33,0.33])
# NB, if you change the above eqn for y, you must regenerate this input:
ax_inset.plot(x_ymax,y_ymax, color = 'm', linestyle = '-', linewidth = '2')
pl.xlabel(r'$\epsilon(0)-1$',size = 14)#\,=\,\frac{2}{\pi}C
pl.ylabel(r'$\rm{MAX}$ $[\zeta/\omega_{0}] $',size = 14)
pl.tick_params(labelsize = 'small')
#ax.legend(loc = 'best')
#pl.savefig('plots/130812_Hopkins_updated_fig_1a.pdf')
pl.savefig('Fig1a.pdf')
pl.savefig('plots/130911_Hopkins_delta_A_max_xi_inset.png', dpi = 300)

pp.savefig()
pl.show()
pp.close()

pl.close()        
