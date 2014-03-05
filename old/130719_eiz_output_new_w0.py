#!/usr/bin/python
import matplotlib               
from numpy import *                    
from pylab import *
from matplotlib import pyplot as pl
from matplotlib import axis as ax
from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('130719_Hopkins_new_peaks_eiz_output.pdf')

## DEFINE FUNCTIONS FOR CALCULATING e(iz)
#---------------------------------------------------------- 
# Matsubara frequencies: z_n at room temp is (2pikbT/hbar)*n (ie coeff*n)
coeff = (2.41*1e14) # in rad/s
n = arange(0,500)
z = n * coeff

# INPUT DATA FILES:
#----------------------------------------------------------
#     WAI YIM: ab initio data for a-sio2 (no peak)
x_nopeak_eV, y_nopeak = loadtxt('a-sio2-eps2.txt', unpack=True, usecols = [0,1])
# DD resent this within set of peak shift data (see x_82_peak_eV):
#x_peak_eV, y_peak = loadtxt('Y24204L.txt', unpack=True, usecols = [0,1])# L => eps2

#----------------------------------------------------------
x_data, y_data = loadtxt('Y24204L.txt', unpack=True, usecols = [0,1])
y_233 =loadtxt('130719_high_w0_y_233.txt', unpack=True, usecols = [0])#peak at 4.2
y_263 =loadtxt('130719_high_w0_y_263.txt', unpack=True, usecols = [0])#peak at 6.2
y_283 =loadtxt('130719_high_w0_y_283.txt', unpack=True, usecols = [0])#peak at 8.2
y_316 =loadtxt('130719_high_w0_y_316.txt', unpack=True, usecols = [0])#peak at 10.2
y_363 =loadtxt('130719_high_w0_y_363.txt', unpack=True, usecols = [0])#peak at 13.2
y_430 =loadtxt('130719_high_w0_y_430.txt', unpack=True, usecols = [0])#peak at 18.2


# MAKE EMPTY LISTS:
#----------------------------------------------------------
x_nopeak  = empty(len(x_nopeak_eV))

eiz_nopeak   = empty(len(z))
eiz_233 = empty(len(z))
eiz_263 = empty(len(z))
eiz_283 = empty(len(z))
eiz_316 = empty(len(z))
eiz_363 = empty(len(z))
eiz_430 = empty(len(z))

diff_eiz_233 =  empty(len(z))
diff_eiz_263 =  empty(len(z))
diff_eiz_283 =  empty(len(z))
diff_eiz_316 =  empty(len(z))
diff_eiz_363 =  empty(len(z))
diff_eiz_430 =  empty(len(z))

# CONVERT ABSCISAE FROM EV TO HERTZ:
#----------------------------------------------------------
x_nopeak  = x_nopeak_eV  *1.5187*10**(15)
x_233  = x_nopeak_eV  *1.5187*10**(15)
x_263  = x_nopeak_eV  *1.5187*10**(15)
x_283  = x_nopeak_eV  *1.5187*10**(15)
x_316  = x_nopeak_eV  *1.5187*10**(15)
x_363  = x_nopeak_eV  *1.5187*10**(15)
x_430  = x_nopeak_eV  *1.5187*10**(15)

## calc difference eps2 and eiz for with and without peak
#-----------------------------------------------------------
diff_eiz_233 = eiz_233 - eiz_nopeak
diff_eiz_263 = eiz_263 - eiz_nopeak
diff_eiz_283 = eiz_283 - eiz_nopeak
diff_eiz_316 = eiz_316 - eiz_nopeak
diff_eiz_363 = eiz_363 - eiz_nopeak
diff_eiz_430 = eiz_430 - eiz_nopeak

listofzeros = np.zeros(len(x_nopeak)) # plot line for y = 0

for j in range(len(z)):
    eiz_nopeak_arg=empty(len(x_nopeak))
    eiz_233_arg=empty(len(x_233))
    eiz_263_arg=empty(len(x_263))
    eiz_283_arg=empty(len(x_283))
    eiz_316_arg=empty(len(x_316))
    eiz_363_arg=empty(len(x_363))
    eiz_430_arg=empty(len(x_430))

    for i in range(len(x_nopeak)):
        eiz_nopeak_arg[i]=x_nopeak[i]*y_nopeak[i] / (x_nopeak[i]**2 + z[j]**2)
    eiz_nopeak[j] = 1 + (2./pi) * trapz(eiz_nopeak_arg,x_nopeak)

    for k in range(len(x_233)):
        eiz_233_arg[k]=x_233[k]*y_233[k]/(x_233[k]**2 + z[j]**2)
    eiz_233[j] = 1 + (2./pi)*trapz(eiz_233_arg,x_233)    

    for m in range(len(x_263)):
        eiz_263_arg[m]=x_263[m]*y_263[m]/(x_263[m]**2 + z[j]**2)
    eiz_263[j] = 1 + (2./pi)*trapz(eiz_263_arg,x_263)    

    for n in range(len(x_283)):
        eiz_283_arg[n]=x_283[n]*y_283[n]/(x_283[n]**2 + z[j]**2)
    eiz_283[j] = 1 + (2./pi)*trapz(eiz_283_arg,x_283)    

    for p in range(len(x_316)):
        eiz_316_arg[p]=x_316[p]*y_316[p]/(x_316[p]**2 + z[j]**2)
    eiz_316[j] = 1 + (2./pi)*trapz(eiz_316_arg,x_316)    

    for q in range(len(x_363)):
        eiz_363_arg[q]=x_363[q]*y_363[q]/(x_363[q]**2 + z[j]**2)
    eiz_363[j] = 1 + (2./pi)*trapz(eiz_363_arg,x_363)    

    for s in range(len(x_430)):
        eiz_430_arg[s]=x_430[s]*y_430[s]/(x_430[s]**2 + z[j]**2)
    eiz_430[j] = 1 + (2./pi)*trapz(eiz_430_arg,x_430)    

savetxt("x_nopeak.txt", x_nopeak)

savetxt("eiz_nopeak.txt", eiz_nopeak)
savetxt("eiz_233.txt", eiz_233)
savetxt("eiz_263.txt", eiz_263)
savetxt("eiz_283.txt", eiz_283)
savetxt("eiz_316.txt", eiz_316)
savetxt("eiz_363.txt", eiz_363)
savetxt("eiz_430.txt", eiz_430)

pl.figure()
pl.plot(x_nopeak, y_nopeak  , label =r'ab initio' )  
pl.plot(x_nopeak, y_233, label =r'$\omega_{233}$ ')
pl.plot(x_nopeak, y_263, label =r'$\omega_{263}$ ')
pl.plot(x_nopeak, y_283, label =r'$\omega_{283}$ ')
pl.plot(x_nopeak, y_316, label =r'$\omega_{316}$ ')
pl.plot(x_nopeak, y_363, label =r'$\omega_{363}$ ')
pl.plot(x_nopeak, y_430, label =r'$\omega_{430}$ ')
pl.xlabel(r'$\omega\,\,[sec^{-1}]$', size = 'x-large')
pl.ylabel(r'$ \epsilon$"($\omega$)', size = 'x-large')
pl.legend(loc = 'best')
pp.savefig()
pl.show()
pl.close

pl.figure()
pl.plot(z, eiz_nopeak, label =r'ab initio')  
pl.plot(z, eiz_233, label =r'$\omega_{233}$')
pl.plot(z, eiz_263, label =r'$\omega_{263}$')
pl.plot(z, eiz_283, label =r'$\omega_{283}$')
pl.plot(z, eiz_316, label =r'$\omega_{316}$')
pl.plot(z, eiz_363, label =r'$\omega_{363}$')
pl.plot(z, eiz_430, label =r'$\omega_{430}$')
pl.xlabel(r'$\xi_{n}$', size = 'x-large')
pl.ylabel(r'$ \epsilon (i \xi_{n} ) $', size = 'x-large')
pl.legend(loc = 'best')
pp.savefig()
pl.show()
pl.close

#pl.figure()
#pl.plot(  z, diff_eiz_233, label =r'$\omega_{042}$')
#pl.plot(  z, diff_eiz_263, label =r'$\omega_{062}$')
#pl.plot(  z, diff_eiz_283, label =r'$\omega_{082}$')
#pl.plot(  z, diff_eiz_316, label =r'$\omega_{102}$')
#pl.plot(  z, diff_eiz_363, label =r'$\omega_{132}$')
#pl.plot(  z, diff_eiz_430, label =r'$\omega_{182}$')
#pl.xlabel(r'$\xi_{n}$', size = 'x-large')
#pl.ylabel(r'$ \delta\epsilon$(i$\xi_{n}$)', size = 'x-large')
#pl.legend(loc = 'best')
#pl.show()
pp.close()










