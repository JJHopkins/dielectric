#!/usr/bin/python
import matplotlib               
from numpy import *                    
from pylab import *
from matplotlib import pyplot as pl
from matplotlib import axis as ax
from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('130715_Hopkins_shifted peaks_eiz.pdf')

## DEFINE FUNCTIONS FOR CALCULATING e(iz)
#---------------------------------------------------------- 
# Matsubara frequencies: z_n at room temp is (2pikbT/hbar)*n (ie coeff*n)
coeff = 0.159 # in eV#(2.41*1e14) # in rad/s
n = arange(0,500)
z = n * coeff

# INPUT DATA FILES:
#----------------------------------------------------------
#     WAI YIM: ab initio data for a-sio2 (no peak)
x_org_nopeak_eV, y_org_nopeak = loadtxt('a-sio2-eps2.txt', unpack=True, usecols = [0,1])
x_nopeak_eV = x_org_nopeak_eV[::3]
y_nopeak = y_org_nopeak[::3]
# DD resent this within set of peak shift data (see x_82_peak_eV):
#x_peak_eV, y_peak = loadtxt('Y24204L.txt', unpack=True, usecols = [0,1])# L => eps2

#     DAN DRYDEN: SYNTH WITH PEAK AT DIFF W_0'S
#----------------------------------------------------------
x_org_182_peak_eV, y_org_182_peak =loadtxt('130711_Y24705L.txt', unpack=True, usecols = [0,1])#peak at 18.2

x_182_peak_eV = x_org_182_peak_eV[::4]
y_182_peak = y_org_182_peak[::4]
# MAKE EMPTY LISTS:
#----------------------------------------------------------
x_nopeak  = empty(len(x_nopeak_eV))
x_182_peak= empty(len(x_182_peak_eV))

eiz_nopeak   = empty(len(z))
eiz_182_peak = empty(len(z))

diff_eiz_182 =  empty(len(z))

# CONVERT ABSCISAE FROM EV TO HERTZ:
#----------------------------------------------------------
x_nopeak  = x_nopeak_eV  #*1.5187*10**(15)
x_182_peak= x_182_peak_eV#*1.5187*10**(15)

listofzeros = np.zeros(len(x_nopeak)) # plot line for y = 0
#pl.figure()
for j in range(len(z)):
    eiz_nopeak_arg=empty(len(x_nopeak))
    eiz_182_peak_arg=empty(len(x_182_peak))
    for i in range(len(x_nopeak)):
        eiz_nopeak_arg[i]=x_nopeak[i]*y_nopeak[i] / (x_nopeak[i]**2 + z[j]**2)
    eiz_nopeak[j] = 1. + (2./pi) * trapz(eiz_nopeak_arg,x_nopeak)

    for s in range(len(x_182_peak)):
        #a = zip(s, x_182_peak[s])
        eiz_182_peak_arg[s]=x_182_peak[s]*y_182_peak[s]/(x_182_peak[s]**2 + z[j]**2)
    eiz_182_peak[j] = 1. + (2./pi)*trapz(eiz_182_peak_arg,x_182_peak)    
    
#pl.plot(x_nopeak,eiz_nopeak_arg)

#pl.plot(x_182_peak,eiz_182_peak_arg)
#pl.close()

pl.figure()
pl.plot(  x_nopeak, y_nopeak  , label =r'ab initio' )  
pl.plot(x_182_peak, y_182_peak, label =r'$\omega_{0} = 2.764$')
#pl.xlabel(r'$\omega\,\,[sec^{-1}]$', size = 'x-large')
#pl.ylabel(r'$\epsilon"(\omega)$', size = 'x-large')
pl.legend(loc = 'best')
pp.savefig()
pl.show()
pl.close


pl.figure()
pl.plot(z, eiz_nopeak, label =r'ab initio')  
pl.plot(z, eiz_182_peak, label =r'$\omega_{0} = 2.764$')
#pl.xlabel(r'$\xi_{n}\,\,[sec^{-1}]$', size = 'x-large')
#pl.ylabel(r'$\epsilon (i \xi_{n} ) $', size = 'x-large')
pl.legend(loc = 'best')
pp.savefig()
pl.show()
pl.close

diff_eiz_182 = eiz_182_peak - eiz_nopeak

pl.figure()
pl.plot(  z, diff_eiz_182, label =r'$\omega_{0} = 2.764$')
#pl.xlabel(r'$\xi_{n}\,\,[sec^{-1}]$', size = 'x-large')
#pl.ylabel(r'$ \delta\epsilon (i \xi_{n} ) $', size = 'x-large')
pl.legend(loc = 'best')
pl.show()
pp.close()
pl.close()











