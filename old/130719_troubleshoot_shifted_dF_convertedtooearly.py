#!/usr/bin/python
import matplotlib               
import numpy as np                  
from pylab import *
from matplotlib import pyplot as pl
from matplotlib import axis as ax
from matplotlib.ticker import MultipleLocator#, FormatStrFormatter
from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('130716_Hopkins_shifted_df_DF.pdf')

x_nopeak = np.loadtxt('x_nopeak.txt', unpack=True, usecols = [0])
y_nopeak = np.loadtxt('a-sio2-eps2.txt', unpack=True, usecols = [1])

x_peak = np.loadtxt('x_peak.txt', unpack=True, usecols = [0])
y_peak = np.loadtxt('Y24204L.txt', unpack=True, usecols = [1])# L => eps2
#
eiz_nopeak = np.loadtxt('eiz_nopeak.txt', unpack=True, usecols = [0])
eiz_peak = np.loadtxt('eiz_peak.txt', unpack=True, usecols = [0])
eiz_042_peak = np.loadtxt('eiz_042_peak.txt', unpack=True, usecols = [0])
eiz_062_peak = np.loadtxt('eiz_062_peak.txt', unpack=True, usecols = [0])
eiz_082_peak = np.loadtxt('eiz_082_peak.txt', unpack=True, usecols = [0])
eiz_102_peak = np.loadtxt('eiz_102_peak.txt', unpack=True, usecols = [0])
eiz_132_peak = np.loadtxt('eiz_132_peak.txt', unpack=True, usecols = [0])
eiz_182_peak = np.loadtxt('eiz_182_peak.txt', unpack=True, usecols = [0])



coeff = (2.41*1e14) # in rad/s
n = arange(0,500)
z = n * coeff
z_y = linspace(0,10,500)
A = 1.0 #[1,4,6]
B = 4.0 #[1,4,6]
C = 6.0 #[1,4,6]

sum_Li2 = np.empty(len(z))
dF_042 = np.empty(len(z))
dF_062 = np.empty(len(z))
dF_082 = np.empty(len(z))
dF_102 = np.empty(len(z))
dF_132 = np.empty(len(z))
dF_182 = np.empty(len(z))
sum_dF_042 = 0.0
sum_dF_062 = 0.0
sum_dF_082 = 0.0
sum_dF_102 = 0.0
sum_dF_132 = 0.0
sum_dF_182 = 0.0

for j in range(len(z)):
    sum_Li2[j] = 0.0
    for m in arange(101) + 1:
        sum_Li2[j] += ((((eiz_nopeak[j] -1)/(eiz_nopeak[j] +1))**2)**m)/(m**2)
    
    dF_042[j] = -sum_Li2[j]\
    *(((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1))**(-1)-((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1)))\
    *(eiz_042_peak[j] - eiz_nopeak[j])/(eiz_nopeak[j])
     
    dF_062[j] = -sum_Li2[j]\
    *(((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1))**(-1)-((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1)))\
    *(eiz_062_peak[j] - eiz_nopeak[j])/(eiz_nopeak[j])
     
    dF_082[j] = -sum_Li2[j]\
    *(((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1))**(-1)-((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1)))\
    *(eiz_082_peak[j] - eiz_nopeak[j])/(eiz_nopeak[j])
     
    dF_102[j] = -sum_Li2[j]\
    *(((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1))**(-1)-((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1)))\
    *(eiz_102_peak[j] - eiz_nopeak[j])/(eiz_nopeak[j])
     
    dF_132[j] = -sum_Li2[j]\
    *(((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1))**(-1)-((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1)))\
    *(eiz_132_peak[j] - eiz_nopeak[j])/(eiz_nopeak[j])
     
    dF_182[j] = -sum_Li2[j]\
    *(((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1))**(-1)-((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1)))\
    *(eiz_182_peak[j] - eiz_nopeak[j])/(eiz_nopeak[j])
     
     
sum_dF_042 = sum(dF_042)
sum_dF_062 = sum(dF_062)
sum_dF_082 = sum(dF_082)
sum_dF_102 = sum(dF_102)
sum_dF_132 = sum(dF_132)
sum_dF_182 = sum(dF_182)

print ('Total dF_042 = %s ' % sum_dF_042) 
print ('Total dF_062 = %s ' % sum_dF_062)
print ('Total dF_082 = %s ' % sum_dF_082)
print ('Total dF_102 = %s ' % sum_dF_102)
print ('Total dF_132 = %s ' % sum_dF_132)
print ('Total dF_182 = %s ' % sum_dF_182)

dF_042[0] = (1./2)*dF_042[0]   
dF_062[0] = (1./2)*dF_062[0]   
dF_082[0] = (1./2)*dF_082[0]   
dF_102[0] = (1./2)*dF_102[0]   
dF_132[0] = (1./2)*dF_132[0]   
dF_182[0] = (1./2)*dF_182[0]   

## calc difference eps2 and eiz for with and without peak
#-----------------------------------------------------------
diff_eps = y_peak - y_nopeak
diff_eiz = eiz_peak - eiz_nopeak
listofzeros = np.zeros(len(x_nopeak)) # plot line for y = 0

## PLOTS
#-------------------------------------------------------------
pl.figure()
pl.plot(n,dF_042,  label = r'$\delta F( \omega_0 = 4.2)$')
pl.plot(n,dF_062,  label = r'$\delta F( \omega_0 = 6.2)$')
pl.plot(n,dF_082,  label = r'$\delta F( \omega_0 = 8.2)$')
pl.plot(n,dF_102,  label = r'$\delta F(\omega_0 = 10.2)$')
pl.plot(n,dF_132,  label = r'$\delta F(\omega_0 = 13.2)$')
pl.plot(n,dF_182,  label = r'$\delta F(\omega_0 = 18.2)$')
#pl.title(r'vdW Free Energy change as a function of n$^{th}$ Matsubara fequency')
pl.tick_params(labelsize = 'small')
pl.xlabel(r'$n$', size = 'large')
pl.ylabel(r'$\delta \mathrm{ F}\,\, \frac{ 8 \pi D^{2} }{k_{B} T}$')
pl.legend(loc = 'best')
pl.text(75,-0.053,r'Total $\delta F( \omega_0 = 4.2) = %.2f$' % sum_dF_042, color='b')
pl.text(75,-0.056,r'Total $\delta F( \omega_0 = 6.2) = %.2f$' % sum_dF_062, color='g')
pl.text(75,-0.059,r'Total $\delta F( \omega_0 = 8.2) = %.2f$' % sum_dF_082, color='r')
pl.text(75,-0.062,r'Total $\delta F(\omega_0 = 10.2) = %.2f$' % sum_dF_102, color='c')
pl.text(75,-0.065,r'Total $\delta F(\omega_0 = 13.2) = %.2f$' % sum_dF_132, color='m')
pl.text(75,-0.068,r'Total $\delta F(\omega_0 = 18.2) = %.2f$' % sum_dF_182, color='y')

pp.savefig()
pl.show()

fig = pl.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8])
ax.plot(x_peak,  y_peak,   color = 'g',label = r'Synthesized')# $\epsilon$"($\omega$)', linestyle='-')
ax.plot(x_nopeak,y_nopeak, color = 'b',label = r'Ab initio')# $\epsilon$"($\omega$)',   linestyle='-')
ax.legend(loc = 'best')
#pl.title(r'$\epsilon$"($\omega$)  Ab Initio and Synthisized')
pl.xlabel(r'$\omega$', size = 'x-large')
pl.ylabel(r'$\epsilon$"($\omega$)', size = 'x-large')

ax_inset = fig.add_axes([0.55,0.42,0.3,0.3])
ax_inset.plot(x_nopeak*1e-16,diff_eps,color='r',label=r'$\epsilon$"($\omega)_{synthesized}$-$\epsilon$"($\omega)_{ab\,initio}$')
ax_inset.plot(x_nopeak*1e-16,listofzeros,color = 'k', label=r'$\delta$$\epsilon$"($\omega$) = 0')
#pl.title(r'Difference $\epsilon$"($\omega$)', size = 'small')
pl.tick_params(labelsize = 'small')
pl.xlabel(r'$\omega$', size = 'small')
pl.ylabel(r'$\delta$$\epsilon$"($\omega$)', size = 'small')
pp.savefig()
pl.show()

fig = pl.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8])
ax.plot(n, eiz_peak, color = 'g', label = r'Synthesized')#$ \epsilon (i \zeta ) _{peak}$')
ax.plot(n, eiz_nopeak, color = 'b', label = r'Ab initio')#$ \epsilon (i \zeta )_{no \,peak} $')
ax.legend(loc = 'best')
#pl.title('$\epsilon$(i$\zeta$) for synthesized and ab initio data')
pl.xlabel(r'$n$', size = 'x-large')
pl.ylabel(r'$ \epsilon (i \zeta ) $', size = 'x-large')

ax_inset = fig.add_axes([0.55,0.42,0.3,0.3])
ax_inset.plot(n,diff_eiz,color= 'r',linestyle='-')
#pl.title('Difference $\epsilon$(i$\zeta$)',size = 'small')
pl.tick_params(labelsize = 'small')
pl.xlabel(r'$n$',size = 'small')#(r'$\zeta$')
pl.ylabel(r'$ \delta \epsilon (i  \zeta )$',size = 'small')
pp.savefig()
pl.show()


pp.close()
pl.close()







