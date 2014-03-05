#!/usr/bin/python
import matplotlib               
import numpy as np                  
from pylab import *
from matplotlib import pyplot as pl
#from matplotlib import axis as ax
from matplotlib.ticker import MultipleLocator#, FormatStrFormatter
from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('plots/130729_Hopkins_shifted_df_DF.pdf')

x_nopeak,y_nopeak=np.loadtxt('data/Y23605L.txt'  , unpack=True, usecols = [0,1])
eiz_nopeak   = np.loadtxt('data/eiz_nopeak.txt'  , unpack=True, usecols = [0])
eiz_042_peak = np.loadtxt('data/eiz_042_peak.txt', unpack=True, usecols = [0])
eiz_062_peak = np.loadtxt('data/eiz_062_peak.txt', unpack=True, usecols = [0])
eiz_082_peak = np.loadtxt('data/eiz_082_peak.txt', unpack=True, usecols = [0])
eiz_102_peak = np.loadtxt('data/eiz_102_peak.txt', unpack=True, usecols = [0])
eiz_132_peak = np.loadtxt('data/eiz_132_peak.txt', unpack=True, usecols = [0])
eiz_182_peak = np.loadtxt('data/eiz_182_peak.txt', unpack=True, usecols = [0])

coeff = 0.159 # in eV. (2.41*1e14) # in rad/s
n = arange(0,250)
z = n * coeff

sum_Li2     = np.empty(len(z))
sum_Li3     = np.empty(len(z))   
sum_Li3_042 = np.empty(len(z))
sum_Li3_062 = np.empty(len(z))
sum_Li3_082 = np.empty(len(z))
sum_Li3_102 = np.empty(len(z))
sum_Li3_132 = np.empty(len(z))
sum_Li3_182 = np.empty(len(z))

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

sum_F     = 0.0
sum_F_042 = 0.0
sum_F_062 = 0.0
sum_F_082 = 0.0
sum_F_102 = 0.0
sum_F_132 = 0.0
sum_F_182 = 0.0

ratio_042 = 0.0
ratio_062 = 0.0
ratio_082 = 0.0
ratio_102 = 0.0
ratio_132 = 0.0
ratio_182 = 0.0

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


    for h in arange(101)+ 1:
        sum_Li3[j] += ((((eiz_nopeak[j] -1)/(eiz_nopeak[j] +1))**2)**h)/(h**3)
        sum_Li3_042[j] += ((((eiz_042_peak[j] -1)/(eiz_042_peak[j] +1))**2)**h)/(h**3)
        sum_Li3_062[j] += ((((eiz_062_peak[j] -1)/(eiz_062_peak[j] +1))**2)**h)/(h**3)
        sum_Li3_082[j] += ((((eiz_082_peak[j] -1)/(eiz_082_peak[j] +1))**2)**h)/(h**3)
        sum_Li3_102[j] += ((((eiz_102_peak[j] -1)/(eiz_102_peak[j] +1))**2)**h)/(h**3)
        sum_Li3_132[j] += ((((eiz_132_peak[j] -1)/(eiz_132_peak[j] +1))**2)**h)/(h**3)
        sum_Li3_182[j] += ((((eiz_182_peak[j] -1)/(eiz_182_peak[j] +1))**2)**h)/(h**3)

dF_042[0] = (1./2)*dF_042[0]   
dF_062[0] = (1./2)*dF_062[0]   
dF_082[0] = (1./2)*dF_082[0]   
dF_102[0] = (1./2)*dF_102[0]   
dF_132[0] = (1./2)*dF_132[0]   
dF_182[0] = (1./2)*dF_182[0]   

sum_Li3[0]= (1./2)*sum_Li3[0]

sum_dF_042 = sum(dF_042)
sum_dF_062 = sum(dF_062)
sum_dF_082 = sum(dF_082)
sum_dF_102 = sum(dF_102)
sum_dF_132 = sum(dF_132)
sum_dF_182 = sum(dF_182)

sum_F     = sum(sum_Li3    )
sum_F_042 = sum(sum_Li3_042)
sum_F_062 = sum(sum_Li3_062)
sum_F_082 = sum(sum_Li3_082)
sum_F_102 = sum(sum_Li3_102)
sum_F_132 = sum(sum_Li3_132)
sum_F_182 = sum(sum_Li3_182)

ratio_042 = 0.0
ratio_062 = 0.0
ratio_082 = 0.0
ratio_102 = 0.0
ratio_132 = 0.0
ratio_182 = 0.0

ratio_042 = sum_dF_042/sum_F
ratio_062 = sum_dF_062/sum_F
ratio_082 = sum_dF_082/sum_F
ratio_102 = sum_dF_102/sum_F
ratio_132 = sum_dF_132/sum_F
ratio_182 = sum_dF_182/sum_F

#diff_sum_F_042 = sum_F- sum_F_042
#diff_sum_F_062 = sum_F- sum_F_062
#diff_sum_F_082 = sum_F- sum_F_082
#diff_sum_F_102 = sum_F- sum_F_102
#diff_sum_F_132 = sum_F- sum_F_132
#diff_sum_F_182 = sum_F- sum_F_182

#perc_diff_df = 0.0
#perc_diff_df = (sum_dF-diff_sum_F)/diff_sum_F

#print ('Total dF   = %s ' % sum_dF)
#print ('Total F_np = %s ' % sum_F)
#print ('Total F_p  = %s ' % sum_F_p)
#print ('Total DF   = %s ' % diff_sum_F)
#print ('Agreement: (dF-DF)/DF = %s ' % perc_diff_df)

#dF[0] = (1./2)*dF[0]   
#sum_Li3[0] = (1./2)*sum_Li3[0]
#sum_Li3_p[0] = (1./2)*sum_Li3_p[0]
#
#F_3 = -sum_Li3
#F_3_p = -sum_Li3_p
#diff_F = (-F_3+F_3_p)
#divide = dF/diff_F 

print ('Total F_nopeak = %s ' % sum_F)

print ('Total dF_042 = %s ' % sum_dF_042) 
print ('Total dF_062 = %s ' % sum_dF_062)
print ('Total dF_082 = %s ' % sum_dF_082)
print ('Total dF_102 = %s ' % sum_dF_102)
print ('Total dF_132 = %s ' % sum_dF_132)
print ('Total dF_182 = %s ' % sum_dF_182)

print ('Ratio dF_042/F_nopeak =  %s ' %ratio_042)
print ('Ratio dF_062/F_nopeak =  %s ' %ratio_062)
print ('Ratio dF_082/F_nopeak =  %s ' %ratio_082)
print ('Ratio dF_102/F_nopeak =  %s ' %ratio_102)
print ('Ratio dF_132/F_nopeak =  %s ' %ratio_132)
print ('Ratio dF_182/F_nopeak =  %s ' %ratio_182)

## calc difference eps2 and eiz for with and without peak
#-----------------------------------------------------------
#diff_eps = y_peak - y_nopeak
#diff_eiz = eiz_peak - eiz_nopeak
#listofzeros = np.zeros(len(x_nopeak)) # plot line for y = 0

## PLOTS
#-------------------------------------------------------------
pl.figure()
pl.plot(n,dF_042, color='g', label = r'$\delta F(\omega_{0} = 0.64)$')
pl.plot(n,dF_062, color='r', label = r'$\delta F(\omega_{0} = 0.94)$')
pl.plot(n,dF_082, color='c', label = r'$\delta F(\omega_{0} = 1.25)$')
pl.plot(n,dF_102, color='m', label = r'$\delta F(\omega_{0} = 1.55)$')
pl.plot(n,dF_132, color='y', label = r'$\delta F(\omega_{0} = 2.01)$')
pl.plot(n,dF_182, color='k', label = r'$\delta F(\omega_{0} = 2.76)$')
#pl.title(r'vdW Free Energy change as a function of n$^{th}$ Matsubara fequency')
pl.tick_params(labelsize = 'small')
pl.xlabel(r'$\xi$', size = 'x-large')
pl.ylabel(r'$\delta \mathrm{ F}\,\, \frac{ 8 \pi D^{2} }{k_{B} T}$', size = 'x-large')
pl.legend(loc = 'best')
#pl.text(75,-0.053,r'$\delta F_{0.64} = %.3f$' % sum_dF_042)#, color='g')
#pl.text(75,-0.056,r'$\delta F_{0.94} = %.3f$' % sum_dF_062)#, color='r')
#pl.text(75,-0.059,r'$\delta F_{1.25} = %.3f$' % sum_dF_082)#, color='c')
#pl.text(75,-0.062,r'$\delta F_{1.55} = %.3f$' % sum_dF_102)#, color='m')
#pl.text(75,-0.065,r'$\delta F_{2.01} = %.3f$' % sum_dF_132)#, color='y')
#pl.text(75,-0.068,r'$\delta F_{2.76} = %.3f$' % sum_dF_182)#, color='k')

pp.savefig()
pl.show()

#fig = pl.figure()
#ax = fig.add_axes([0.1,0.1,0.8,0.8])
##ax.plot(x_peak,  y_peak,   color = 'g',label = r'Synthesized')# $\epsilon$"($\omega$)', linestyle='-')
#ax.plot(x_nopeak,y_nopeak, color = 'b',label = r'Ab initio')# $\epsilon$"($\omega$)',   linestyle='-')
#ax.legend(loc = 'best')
##pl.title(r'$\epsilon$"($\omega$)  Ab Initio and Synthisized')
#pl.xlabel(r'$\omega$', size = 'x-large')
#pl.ylabel(r'$\epsilon$"($\omega$)', size = 'x-large')
#
#ax_inset = fig.add_axes([0.55,0.42,0.3,0.3])
##ax_inset.plot(x_nopeak,diff_eps,color='r',label=r'$\epsilon$"($\omega)_{synthesized}$-$\epsilon$"($\omega)_{ab\,initio}$')
##ax_inset.plot(x_nopeak*1e-16,diff_eps,color='r',label=r'$\epsilon$"($\omega)_{synthesized}$-$\epsilon$"($\omega)_{ab\,initio}$')
##ax_inset.plot(x_nopeak*1e-16,listofzeros,color = 'k', label=r'$\delta$$\epsilon$"($\omega$) = 0')
##ax_inset.plot(x_nopeak,listofzeros,color = 'k', label=r'$\delta$$\epsilon$"($\omega$) = 0')
##pl.title(r'Difference $\epsilon$"($\omega$)', size = 'small')
#pl.tick_params(labelsize = 'small')
#pl.xlabel(r'$\omega$', size = 'small')
#pl.ylabel(r'$\delta$$\epsilon$"($\omega$)', size = 'small')
#pp.savefig()
#pl.show()

fig = pl.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8])
ax.plot(z,  eiz_nopeak, color='b', label = r'$Ab\, initio$')#$ \epsilon (i \zeta )_{no \,peak} $')
ax.plot(z,eiz_042_peak, color='g', label = r'$\omega_{0} = 0.64$')
ax.plot(z,eiz_062_peak, color='r', label = r'$\omega_{0} = 0.94$')
ax.plot(z,eiz_082_peak, color='c', label = r'$\omega_{0} = 1.25$')
ax.plot(z,eiz_102_peak, color='m', label = r'$\omega_{0} = 1.55$')
ax.plot(z,eiz_132_peak, color='y', label = r'$\omega_{0} = 2.01$')
ax.plot(z,eiz_182_peak, color='k', label = r'$\omega_{0} = 2.76$')

ax.legend(loc = 'lower right')
#pl.title('$\epsilon$(i$\zeta$) for synthesized and ab initio data')
#pl.xlabel(r'$\frac{\xi_{n}}{\omega_{0}}$', size = 'x-large')
pl.xlabel(r'$\xi$', size = 'x-large')
pl.ylabel(r'$ \epsilon (i \zeta ) $', size = 'x-large')

diff_eiz_042 = eiz_042_peak - eiz_nopeak
diff_eiz_062 = eiz_062_peak - eiz_nopeak
diff_eiz_082 = eiz_082_peak - eiz_nopeak
diff_eiz_102 = eiz_102_peak - eiz_nopeak
diff_eiz_132 = eiz_132_peak - eiz_nopeak
diff_eiz_182 = eiz_182_peak - eiz_nopeak

ax_inset = fig.add_axes([0.45,0.55,0.33,0.33])
ax_inset.plot(z,diff_eiz_042,color= 'g')#,linestyle='--')
ax_inset.plot(z,diff_eiz_062,color= 'r')#,linestyle='--')
ax_inset.plot(z,diff_eiz_082,color= 'c')#,linestyle='--')
ax_inset.plot(z,diff_eiz_102,color= 'm')#,linestyle='--')
ax_inset.plot(z,diff_eiz_132,color= 'y')#,linestyle='--')
ax_inset.plot(z,diff_eiz_182,color= 'k')#,linestyle='--')
#pl.title('Difference $\epsilon$(i$\zeta$)',size = 'small')
pl.tick_params(labelsize = 'small')
pl.xlabel(r'$\xi$',size = 'small')#(r'$\zeta$')
#pl.xlabel(r'$\frac{\xi_{n}}{\omega_{0}}$', size = 'small')
pl.ylabel(r'$ \delta \epsilon (i  \zeta )$',size = 'small')
pp.savefig()
pl.show()

pl.close()
pp.close()









