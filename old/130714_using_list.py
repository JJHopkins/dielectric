#!/usr/bin/python
import matplotlib               
import numpy                    
from pylab import *
from matplotlib import pyplot as pl
from matplotlib import axis as ax
from matplotlib.ticker import MultipleLocator#, FormatStrFormatter
#from scipy import integrate
from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('130711_Hopkins_halved_df_DF.pdf')

x_nopeak = numpy.loadtxt('x_nopeak.txt', unpack=True, usecols = [0])
y_nopeak = numpy.loadtxt('a-sio2-eps2.txt', unpack=True, usecols = [1])

x_peak = numpy.loadtxt('x_peak.txt', unpack=True, usecols = [0])
y_peak = numpy.loadtxt('Y24204L.txt', unpack=True, usecols = [1])# L => eps2

x_Re_expt_eV, y_Re_expt = numpy.loadtxt('Y24204K.txt', unpack=True, usecols = [0,1])# K => eps1


eiz_nopeak = numpy.loadtxt('eiz_nopeak.txt', unpack=True, usecols = [0])
eiz_peak = numpy.loadtxt('eiz_peak.txt', unpack=True, usecols = [0])

coeff = (2.41*1e14) # in rad/s
n = arange(0,250)
z = n * coeff
w_0 = 1e16

z_ys = linspace(0,2e16,250)
As = [1,2,4,6,8,10,12]
sum_Li2_y = numpy.empty(len(z_ys))
qs = range(1,101)
def Lamda(a,xi):
    return a * xi**2 * (xi/w_0) /\
    (2* xi**3 (1 + (xi/w_0)**2) + a*(xi/w_0)**3)

def Li2(lamda,power):
    return (lamda**2)**power/\
    (power**2) 

def y(dilog, lamda):
    return dilog*\
    (1 - lamda**2)/\
    lamda

for j, A in enumerate(As):
    for i, z_y in enumerate(z_ys):
	contrast = Lamda(A, z_y)
        for k in enumerate(qs):
	    sum_Li2_y += [Li2(contrast, q) for q in qs]
        ys = [y(sum_Li2_y, contrast) for A in As]
        ys_normed = [ys/y[0] for A in As]
	max_ys = [max(ys_normed) for A in As] 

sum_Li2 = numpy.empty(len(z))
dF = numpy.empty(len(z))
sum_dF =  0.0
sum_F =  0.0
sum_F_p =  0.0
diff_sum_F =  0.0

sum_Li3 = numpy.empty(len(z))
sum_Li3_p = numpy.empty(len(z))
F_3 = numpy.empty(len(z))
F_3_p = numpy.empty(len(z))

pl.figure()
pl.plot(z_ys/w_0,ys_normed)
pl.title(r'$Function\,\,y(\xi)\,\,from \,\,eqn\,6$')
pl.tick_params(labelsize = 'small')
pl.xlabel(r'$0<\frac{\xi}{\omega_0}<2$', size = 'large')
pl.ylabel(r'$y(\xi)$')# \mathrm{ F}\,\, \frac{ 8 \pi D^{2} }{k_{B} T}$')
pl.legend(loc = 'best')
pl.show()

pl.figure()
pl.plot(As,max_ys)
pl.title(r'max values for y(xi)')
pl.tick_params(labelsize = 'small')
pl.xlabel(r'As')#$0<\frac{\xi}{\omega_0}<2$', size = 'large')
pl.ylabel(r'max y')#$y(\xi)$')# \mathrm{ F}\,\, \frac{ 8 \pi D^{2} }{k_{B} T}$')
pl.legend(loc = 'best')
pl.show()

for j in range(len(z)):
    sum_Li2[j] = 0.0
    for m in arange(101) + 1:
        sum_Li2[j] += ((((eiz_nopeak[j] -1)/(eiz_nopeak[j] +1))**2)**m)/(m**2)
    dF[j] = -sum_Li2[j]\
    *(((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1))**(-1)-((eiz_nopeak[j] -1)/(eiz_nopeak[j] + 1)))\
    *(eiz_peak[j] - eiz_nopeak[j])/(eiz_nopeak[j])
     
    sum_Li3[j] = 0.0
    sum_Li3_p[j] = 0.0
    
    for h in arange(101)+ 1:
        sum_Li3[j] += ((((eiz_nopeak[j] -1)/(eiz_nopeak[j] +1))**2)**h)/(h**3)
        sum_Li3_p[j] += ((((eiz_peak[j] -1)/(eiz_peak[j] +1))**2)**h)/(h**3)


sum_dF = sum(dF)
sum_F = sum(sum_Li3)
sum_F_p = sum(sum_Li3_p)
diff_sum_F = sum_F- sum_F_p
perc_diff_df = 0.0
perc_diff_df = (sum_dF-diff_sum_F)/diff_sum_F

print ('Total dF   = %s ' % sum_dF)
print ('Total F_np = %s ' % sum_F)
print ('Total F_p  = %s ' % sum_F_p)
print ('Total DF   = %s ' % diff_sum_F)
print ('Agreement: (dF-DF)/DF = %s ' % perc_diff_df)

dF[0] = (1./2)*dF[0]   
sum_Li3[0] = (1./2)*sum_Li3[0]
sum_Li3_p[0] = (1./2)*sum_Li3_p[0]

F_3 = -sum_Li3
F_3_p = -sum_Li3_p
diff_F = (-F_3+F_3_p)
divide = dF/diff_F 

## calc difference eps2 and eiz for with and without peak
#-----------------------------------------------------------
diff_eps = y_peak - y_nopeak
diff_eiz = eiz_peak - eiz_nopeak
listofzeros = numpy.zeros(len(x_nopeak)) # plot line for y = 0

## PLOTS
#-------------------------------------------------------------
fig = pl.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8])
ax.plot(x_peak,  y_peak,   color = 'g',label = r'Synthesized')# $\epsilon$"($\omega$)', linestyle='-')
ax.plot(x_nopeak,y_nopeak, color = 'b',label = r'Ab initio')# $\epsilon$"($\omega$)',   linestyle='-')
ax.legend(loc = 'best')
#pl.title(r'$\epsilon$"($\omega$)  Ab Initio and Synthisized')
pl.xlabel(r'$\omega$', size = 'large')
pl.ylabel(r'$\epsilon$"($\omega$)')

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
pl.xlabel(r'$n$', size = 'large')
pl.ylabel(r'$ \epsilon (i \zeta ) $')

ax_inset = fig.add_axes([0.55,0.42,0.3,0.3])
ax_inset.plot(n,diff_eiz,color= 'r',linestyle='-')
#pl.title('Difference $\epsilon$(i$\zeta$)',size = 'small')
pl.tick_params(labelsize = 'small')
pl.xlabel(r'$n$',size = 'small')#(r'$\zeta$')
pl.ylabel(r'$ \delta \epsilon (i  \zeta )$',size = 'small')
pp.savefig()
pl.show()

pl.figure()
pl.plot(n,F_3, color = 'k', label = r'$\left( F/S \right)_{ab\,initio}$')
pl.plot(n,F_3_p,color='k',linestyle='--',label=r'$\left( F/S \right) _{synthesized}$')
#pl.title(r'vdW Free Energy as a function of n$^{th}$ Matsubara fequency')
pl.tick_params(labelsize = 'small')
pl.xlabel(r'$n$', size = 'large')
pl.ylabel(r'$\left( \mathrm{F/S} \right) \,\, \frac{ 8 \pi D^{2} } {k _{B} T}$')
pl.legend(loc='best')
pp.savefig()
pl.show()

pl.figure()
pl.plot(n,dF, color = 'r', label = r'$\delta F$')
pl.plot(n,diff_F, color = 'k', linestyle = '--', label = r'$\Delta F$')
#pl.title(r'vdW Free Energy change as a function of n$^{th}$ Matsubara fequency')
pl.tick_params(labelsize = 'small')
pl.xlabel(r'$n$', size = 'large')
pl.ylabel(r'$\delta \mathrm{ F}\,\, \frac{ 8 \pi D^{2} }{k_{B} T}$')
pl.legend(loc = 'best')
pp.savefig()
pl.show()

pp.close()
pl.close()





