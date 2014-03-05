#!/usr/bin/python
import matplotlib               
import numpy as np                  
from pylab import *
from matplotlib import pyplot as pl
from matplotlib import axis as ax
from matplotlib.ticker import MultipleLocator#, FormatStrFormatter
from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('plots/130729_Hopkins_new_w0_df_DF.pdf')

x_nopeak = np.loadtxt('data/x_nopeak.txt', unpack=True, usecols = [0])
y_nopeak = np.loadtxt('data/a-sio2-eps2.txt', unpack=True, usecols = [1])

x_peak = np.loadtxt('data/x_peak.txt', unpack=True, usecols = [0])
y_peak = np.loadtxt('data/Y24204L.txt', unpack=True, usecols = [1])# L => eps2
#
eiz_nopeak = np.loadtxt('data/eiz_nopeak.txt', unpack=True, usecols = [0])
eiz_peak = np.loadtxt('data/eiz_peak.txt', unpack=True, usecols = [0])
eiz_233 = np.loadtxt('data/eiz_233.txt', unpack=True, usecols = [0])
eiz_263 = np.loadtxt('data/eiz_263.txt', unpack=True, usecols = [0])
eiz_283 = np.loadtxt('data/eiz_283.txt', unpack=True, usecols = [0])
eiz_316 = np.loadtxt('data/eiz_316.txt', unpack=True, usecols = [0])
eiz_363 = np.loadtxt('data/eiz_363.txt', unpack=True, usecols = [0])
eiz_430 = np.loadtxt('data/eiz_430.txt', unpack=True, usecols = [0])



coeff = (2.41*1e14) # in rad/s
n = arange(0,500)
z = n * coeff
z_y = linspace(0,10,500)
A = 1.0 #[1,4,6]
B = 4.0 #[1,4,6]
C = 6.0 #[1,4,6]

sum_Li2 = np.empty(len(z))
dF_233= np.empty(len(z))
dF_263= np.empty(len(z))
dF_283= np.empty(len(z))
dF_316= np.empty(len(z))
dF_363= np.empty(len(z))
dF_430= np.empty(len(z))
sum_dF_233 = 0.0
sum_dF_263 = 0.0
sum_dF_283 = 0.0
sum_dF_316 = 0.0
sum_dF_363 = 0.0
sum_dF_430 = 0.0

for j in range(len(z)):
    sum_Li2[j] = 0.0
    for m in arange(101) + 1:
        sum_Li2[j] += ((((eiz_nopeak[j] -1)/(eiz_nopeak[j] +1))**2)**m)/(m**2)
    
    dF_233[j] = -sum_Li2[j]\
    *(((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1))**(-1)-((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1)))\
    *(eiz_233[j] - eiz_nopeak[j])/(eiz_nopeak[j])
     
    dF_263[j] = -sum_Li2[j]\
    *(((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1))**(-1)-((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1)))\
    *(eiz_263[j] - eiz_nopeak[j])/(eiz_nopeak[j])
     
    dF_283[j] = -sum_Li2[j]\
    *(((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1))**(-1)-((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1)))\
    *(eiz_283[j] - eiz_nopeak[j])/(eiz_nopeak[j])
     
    dF_316[j] = -sum_Li2[j]\
    *(((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1))**(-1)-((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1)))\
    *(eiz_316[j] - eiz_nopeak[j])/(eiz_nopeak[j])
     
    dF_363[j] = -sum_Li2[j]\
    *(((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1))**(-1)-((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1)))\
    *(eiz_363[j] - eiz_nopeak[j])/(eiz_nopeak[j])
     
    dF_430[j] = -sum_Li2[j]\
    *(((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1))**(-1)-((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1)))\
    *(eiz_430[j] - eiz_nopeak[j])/(eiz_nopeak[j])
     
     
sum_dF_233 = sum(dF_233)
sum_dF_263 = sum(dF_263)
sum_dF_283 = sum(dF_283)
sum_dF_316 = sum(dF_316)
sum_dF_363 = sum(dF_363)
sum_dF_430 = sum(dF_430)

print ('Total dF_233 = %s ' % sum_dF_233) 
print ('Total dF_263 = %s ' % sum_dF_263)
print ('Total dF_283 = %s ' % sum_dF_283)
print ('Total dF_316 = %s ' % sum_dF_316)
print ('Total dF_363 = %s ' % sum_dF_363)
print ('Total dF_430 = %s ' % sum_dF_430)

dF_233[0] = (1./2)*dF_233[0]   
dF_263[0] = (1./2)*dF_263[0]   
dF_283[0] = (1./2)*dF_283[0]   
dF_316[0] = (1./2)*dF_316[0]   
dF_363[0] = (1./2)*dF_363[0]   
dF_430[0] = (1./2)*dF_430[0]   

## calc difference eps2 and eiz for with and without peak
#-----------------------------------------------------------
diff_eps = y_peak - y_nopeak
diff_eiz = eiz_peak - eiz_nopeak
listofzeros = np.zeros(len(x_nopeak)) # plot line for y = 0

## PLOTS
#-------------------------------------------------------------
pl.figure()
pl.plot(n,dF_233, color='g',  label = r'$\delta F(\omega_{233})$')
pl.plot(n,dF_263, color='r',  label = r'$\delta F(\omega_{263})$')
pl.plot(n,dF_283, color='c',  label = r'$\delta F(\omega_{283})$')
pl.plot(n,dF_316, color='m',  label = r'$\delta F(\omega_{316})$')
pl.plot(n,dF_363, color='y',  label = r'$\delta F(\omega_{363})$')
pl.plot(n,dF_430, color='k',  label = r'$\delta F(\omega_{430})$')
#pl.title(r'vdW Free Energy change as a function of n$^{th}$ Matsubara fequency')
pl.tick_params(labelsize = 'small')
pl.xlabel(r'$n$', size = 'large')
pl.ylabel(r'$\delta \mathrm{ F}\,\, \frac{ 8 \pi D^{2} }{k_{B} T}$')
pl.legend(loc = 'best')
pl.text(200,-0.0040,r'Total $\delta F(\omega_{233}) = %.2f$' % sum_dF_233)#, color='g')
pl.text(200,-0.0044,r'Total $\delta F(\omega_{263}) = %.2f$' % sum_dF_263)#, color='r')
pl.text(200,-0.0048,r'Total $\delta F(\omega_{283}) = %.2f$' % sum_dF_283)#, color='c')
pl.text(200,-0.0052,r'Total $\delta F(\omega_{316}) = %.2f$' % sum_dF_316)#, color='m')
pl.text(200,-0.0056,r'Total $\delta F(\omega_{363}) = %.2f$' % sum_dF_363)#, color='y')
pl.text(200,-0.0060,r'Total $\delta F(\omega_{430}) = %.2f$' % sum_dF_430)#, color='k')

pp.savefig()
pl.show()

#fig = pl.figure()
#ax = fig.add_axes([0.1,0.1,0.8,0.8])
#ax.plot(x_peak,  y_peak,  color= 'g',label=r'Synthesized')# $\epsilon$"($\omega$)',linestyle='-')
#ax.plot(x_nopeak,y_nopeak,color= 'b',label=r'Ab initio')# $\epsilon$"($\omega$)',linestyle='-')
#ax.legend(loc = 'best')
##pl.title(r'$\epsilon$"($\omega$)  Ab Initio and Synthisized')
#pl.xlabel(r'$\omega$', size = 'x-large')
#pl.ylabel(r'$\epsilon$"($\omega$)', size = 'x-large')
#
#ax_inset = fig.add_axes([0.55,0.42,0.3,0.3])
#ax_inset.plot(x_nopeak*1e-16,diff_eps,color='r',label=r'$\epsilon$"($\omega)_{synthesized}$-$\epsilon$"($\omega)_{ab\,initio}$')
#ax_inset.plot(x_nopeak*1e-16,listofzeros,color = 'k', label=r'$\delta$$\epsilon$"($\omega$) = 0')
##pl.title(r'Difference $\epsilon$"($\omega$)', size = 'small')
#pl.tick_params(labelsize = 'small')
#pl.xlabel(r'$\omega$', size = 'small')
#pl.ylabel(r'$\delta$$\epsilon$"($\omega$)', size = 'small')
#pp.savefig()
#pl.show()
#
#fig = pl.figure()
#ax = fig.add_axes([0.1,0.1,0.8,0.8])
#ax.plot(n, eiz_peak, color = 'g', label = r'Synthesized')#$ \epsilon (i \zeta ) _{peak}$')
#ax.plot(n, eiz_nopeak, color = 'b', label = r'Ab initio')#$ \epsilon (i \zeta )_{no \,peak} $')
#ax.legend(loc = 'best')
##pl.title('$\epsilon$(i$\zeta$) for synthesized and ab initio data')
#pl.xlabel(r'$n$', size = 'x-large')
#pl.ylabel(r'$ \epsilon (i \zeta ) $', size = 'x-large')
#
#ax_inset = fig.add_axes([0.55,0.42,0.3,0.3])
#ax_inset.plot(n,diff_eiz,color= 'r',linestyle='-')
##pl.title('Difference $\epsilon$(i$\zeta$)',size = 'small')
#pl.tick_params(labelsize = 'small')
#pl.xlabel(r'$n$',size = 'small')#(r'$\zeta$')
#pl.ylabel(r'$ \delta \epsilon (i  \zeta )$',size = 'small')
#pp.savefig()
#pl.show()


pp.close()
pl.close()







