#!/usr/bin/python
import matplotlib               
from numpy import *                    
from pylab import *
from matplotlib import pyplot as pl
from matplotlib import axis as ax
from matplotlib.backends.backend_pdf import PdfPages
#pp = PdfPages('plots/130729_Hopkins_shifted peaks_eiz.pdf')

## DEFINE FUNCTIONS FOR CALCULATING e(iz)
#---------------------------------------------------------- 
# Matsubara frequencies: z_n at room temp is (2pikbT/hbar)*n (ie coeff*n)
coeff = 0.159 # in eV#(2.41*1e14) # in rad/s
n = arange(0,500)
z = n * coeff

# INPUT DATA FILES:
#----------------------------------------------------------
#     WAI YIM: ab initio data for a-sio2 (no peak)
x_nopeak_eV, y_nopeak = loadtxt('data/Y23605L.txt', unpack=True, usecols = [0,1])
# DD resent this within set of peak shift data (see x_82_peak_eV):
#x_peak_eV, y_peak = loadtxt('Y24204L.txt', unpack=True, usecols = [0,1])# L => eps2

#     DAN DRYDEN: SYNTH WITH PEAK AT DIFF W_0'S
#----------------------------------------------------------
x_042_peak_eV,y_042_peak_eV=loadtxt('data/CollEt.csv', delimiter = ',',unpack=True, usecols = [0,1])#peak at 4.2
x_062_peak_eV,y_062_peak_eV=loadtxt('data/CollEx.csv', delimiter = ',',unpack=True, usecols = [0,1])#peak at 6.2
x_082_peak_eV,y_082_peak_eV=loadtxt('data/CollEy.csv', delimiter = ',',unpack=True, usecols = [0,1])#peak at 8.2
x_102_peak_eV,y_102_peak_eV=loadtxt('data/CollEz.csv', delimiter = ',',unpack=True, usecols= [0,1])#peak at 10.2

# MAKE EMPTY LISTS:
#----------------------------------------------------------
x_042_peak= empty(len(x_042_peak_eV))
x_062_peak= empty(len(x_062_peak_eV))
x_082_peak= empty(len(x_082_peak_eV))
x_102_peak= empty(len(x_102_peak_eV))

y_042_peak= empty(len(y_042_peak_eV))
y_062_peak= empty(len(y_062_peak_eV))
y_082_peak= empty(len(y_082_peak_eV))
y_102_peak= empty(len(y_102_peak_eV))

eiz_042_peak = empty(len(z))
eiz_062_peak = empty(len(z))
eiz_082_peak = empty(len(z))
eiz_102_peak = empty(len(z))

# CONVERT ABSCISAE FROM EV TO HERTZ:
#----------------------------------------------------------
#x_nopeak  = x_nopeak_eV  #*1.5187*10**(15)
x_042_peak= x_042_peak_eV#*1.5187*10**(15)
x_062_peak= x_062_peak_eV#*1.5187*10**(15)
x_082_peak= x_082_peak_eV#*1.5187*10**(15)
x_102_peak= x_102_peak_eV#*1.5187*10**(15)

y_042_peak= 9.3929*y_042_peak_eV#*1.5187*10**(15)
y_062_peak= 9.3929*y_062_peak_eV#*1.5187*10**(15)
y_082_peak= 9.3929*y_082_peak_eV#*1.5187*10**(15)
y_102_peak= 9.3929*y_102_peak_eV#*1.5187*10**(15)

listofzeros = np.zeros(len(x_042_peak)) # plot line for y = 0
for j in range(len(z)):
    eiz_042_peak_arg=empty(len(x_042_peak))
    eiz_062_peak_arg=empty(len(x_062_peak))
    eiz_082_peak_arg=empty(len(x_082_peak))
    eiz_102_peak_arg=empty(len(x_102_peak))

    for k in range(len(x_042_peak)):
        eiz_042_peak_arg[k]=x_042_peak[k]*y_042_peak[k]/(x_042_peak[k]**2 + z[j]**2)
    eiz_042_peak[j] = 1. + (2./pi)*trapz(eiz_042_peak_arg,x_042_peak)    

    for m in range(len(x_062_peak)):
        eiz_062_peak_arg[m]=x_062_peak[m]*y_062_peak[m]/(x_062_peak[m]**2 + z[j]**2)
    eiz_062_peak[j] = 1. + (2./pi)*trapz(eiz_062_peak_arg,x_062_peak)    

    for n in range(len(x_082_peak)):
        eiz_082_peak_arg[n]=x_082_peak[n]*y_082_peak[n]/(x_082_peak[n]**2 + z[j]**2)
    eiz_082_peak[j] = 1. + (2./pi)*trapz(eiz_082_peak_arg,x_082_peak)    

    for p in range(len(x_102_peak)):
        eiz_102_peak_arg[p]=x_102_peak[p]*y_102_peak[p]/(x_102_peak[p]**2 + z[j]**2)
    eiz_102_peak[j] = 1. + (2./pi)*trapz(eiz_102_peak_arg,x_102_peak)    

pl.figure()
pl.plot(x_042_peak, y_042_peak,color = 'k', label =r'$\epsilon$"$(\omega)$')
pl.plot(x_062_peak, y_062_peak,color = 'b', label =r'$\epsilon_{\hat{x}}$"$(\omega)$')
pl.plot(x_082_peak, y_082_peak,color = 'g', label =r'$\epsilon_{\hat{y}}$"$(\omega)$')
pl.plot(x_102_peak, y_102_peak,color = 'r', label =r'$\epsilon_{\hat{z}}$"$(\omega)$')
pl.xlabel(r'$\hbar\omega\,\,[eV]$', size = 21)
#pl.legend(loc = 'best')
pl.ylabel(r'$\epsilon$"($\omega$)', size = 21)
pl.legend(loc = 'best')
#pl.axis([0,35,0,0.27])
#pl.savefig('plots/130814_collagen_eps2.pdf')
#pl.savefig('plots/130814_collagen_eps2.jpg', dpi = 300)
pl.savefig('plots/130815_collagen_eps2.png', dpi = 300)

pl.figure()
#pl.plot(z, eiz_nopeak, label =r'$ab\, initio$')  
pl.plot(z, eiz_042_peak,color = 'k')#, label =r'$\omega_{0} = 0.64$') 
pl.plot(z, eiz_062_peak,color = 'b')#, label =r'$\omega_{0} = 0.94$') 
pl.plot(z, eiz_082_peak,color = 'g')#, label =r'$\omega_{0} = 1.25$') 
pl.plot(z, eiz_102_peak,color = 'r')#, label =r'$\omega_{0} = 1.55$')
pl.xlabel(r'$\xi_{N}\,\,[eV]$', size = 21)#'x-large')
pl.ylabel(r'$\epsilon (i \xi_{N} ) $', size = 21)#'x-large')
pl.axis([0,35,1.15,2.5])
#pl.title(r'$\epsilon(i\xi)\, \rm{for\, ab\, initio\, and\, L.O.\, peaks \,at}\, \omega_{0}$')
#pl.legend(loc = 'best')
#pl.savefig('plots/130814_collagen_eiz.pdf')
#pl.savefig('plots/130814_collagen_eiz.jpg', dpi = 300)
pl.savefig('plots/130815_collagen_eiz.png', dpi = 300)
pl.show()
pl.close()


#*10**(-16)
#*10**(-16)
#*10**(-16)
#*10**(-16)




