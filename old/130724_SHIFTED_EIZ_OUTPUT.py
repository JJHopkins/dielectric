#!/usr/bin/python
import matplotlib               
from numpy import *                    
from pylab import *
from matplotlib import pyplot as pl
from matplotlib import axis as ax
from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('plots/130724_Hopkins_shifted peaks_eiz.pdf')

## DEFINE FUNCTIONS FOR CALCULATING e(iz)
#---------------------------------------------------------- 
# Matsubara frequencies: z_n at room temp is (2pikbT/hbar)*n (ie coeff*n)
coeff = 0.159 # in eV#(2.41*1e14) # in rad/s
n = arange(0,500)
z = n * coeff

# INPUT DATA FILES:
#----------------------------------------------------------
#     WAI YIM: ab initio data for a-sio2 (no peak)
x_nopeak_eV, y_nopeak = loadtxt('130723_data/Y23605L.txt', unpack=True, usecols = [0,1])
# DD resent this within set of peak shift data (see x_82_peak_eV):
#x_peak_eV, y_peak = loadtxt('Y24204L.txt', unpack=True, usecols = [0,1])# L => eps2

#     DAN DRYDEN: SYNTH WITH PEAK AT DIFF W_0'S
#----------------------------------------------------------
x_042_peak_eV,y_042_peak=loadtxt('130723_data/Y24901L.txt',unpack=True, usecols = [0,1])#peak at 4.2
x_062_peak_eV,y_062_peak=loadtxt('130723_data/Y24902L.txt',unpack=True, usecols = [0,1])#peak at 6.2
x_082_peak_eV,y_082_peak=loadtxt('130723_data/Y24903L.txt',unpack=True, usecols = [0,1])#peak at 8.2
x_102_peak_eV,y_102_peak=loadtxt('130723_data/Y24904L.txt',unpack=True, usecols= [0,1])#peak at 10.2
x_132_peak_eV,y_132_peak=loadtxt('130723_data/Y24905L.txt',unpack=True, usecols= [0,1])#peak at 13.2
x_182_peak_eV,y_182_peak=loadtxt('130723_data/Y24906L.txt',unpack=True, usecols= [0,1])#peak at 18.2


x_s,y_s=loadtxt('130723_data/Y24204S.PRN',unpack=True, usecols = [0,1])#peak at 4.2
x_042_s,y_042_s=loadtxt('130723_data/Y24901S.PRN',unpack=True, usecols = [0,1])#peak at 4.2
x_062_s,y_062_s=loadtxt('130723_data/Y24902S.PRN',unpack=True, usecols = [0,1])#peak at 6.2
x_082_s,y_082_s=loadtxt('130723_data/Y24903S.PRN',unpack=True, usecols = [0,1])#peak at 8.2
x_102_s,y_102_s=loadtxt('130723_data/Y24904S.PRN',unpack=True, usecols= [0,1])#peak at 10.2
x_132_s,y_132_s=loadtxt('130723_data/Y24905S.PRN',unpack=True, usecols= [0,1])#peak at 13.2
x_182_s,y_182_s=loadtxt('130723_data/Y24906S.PRN',unpack=True, usecols= [0,1])#peak at 18.2

pl.figure()
pl.plot(x_s,y_s,  label = r'$ab\,initio$')
pl.plot(x_042_s,y_042_s,  label = r'$\omega_0 = 4.2$')
pl.plot(x_062_s,y_062_s,  label = r'$\omega_0 = 6.2$')
pl.plot(x_082_s,y_082_s,  label = r'$\omega_0 = 8.2$')
pl.plot(x_102_s,y_102_s,  label = r'$\omega_0 = 10.2$')
pl.plot(x_132_s,y_132_s,  label = r'$\omega_0 = 13.2$')
pl.plot(x_182_s,y_182_s,  label = r'$\omega_0 = 18.2$')
pl.xlabel(r'$eV$')#, size = 'x-large')
pl.ylabel(r'$Strength$')#, size = 'x-large')
pl.title('Sum Rule')
pl.legend(loc='lower right')
pl.show()
pl.close()

I_ab = trapz(y_s,x_s)
# 285.14086243166378
I_042 = trapz(y_042_s,x_042_s)
# 286.19703515281657
I_062 = trapz(y_062_s,x_062_s)
# 285.67258664466283
I_082 = trapz(y_082_s,x_082_s)
# 285.12940733546276
I_102 = trapz(y_102_s,x_102_s)
# 284.57648538356881
I_132 = trapz(y_132_s,x_132_s)
# 283.75073833021582
I_182 = trapz(y_182_s,x_182_s)
# 282.36586632109379
print (I_042-I_ab)/I_ab
print (I_062-I_ab)/I_ab
print (I_082-I_ab)/I_ab
print (I_102-I_ab)/I_ab
print (I_132-I_ab)/I_ab
print (I_182-I_ab)/I_ab
# MAKE EMPTY LISTS:
#----------------------------------------------------------
x_nopeak  = empty(len(x_nopeak_eV))
x_042_peak= empty(len(x_042_peak_eV))
x_062_peak= empty(len(x_062_peak_eV))
x_082_peak= empty(len(x_082_peak_eV))
x_102_peak= empty(len(x_102_peak_eV))
x_132_peak= empty(len(x_132_peak_eV))
x_182_peak= empty(len(x_182_peak_eV))

eiz_nopeak   = empty(len(z))
eiz_042_peak = empty(len(z))
eiz_062_peak = empty(len(z))
eiz_082_peak = empty(len(z))
eiz_102_peak = empty(len(z))
eiz_132_peak = empty(len(z))
eiz_182_peak = empty(len(z))

diff_eiz_042 =  empty(len(z))
diff_eiz_062 =  empty(len(z))
diff_eiz_082 =  empty(len(z))
diff_eiz_102 =  empty(len(z))
diff_eiz_132 =  empty(len(z))
diff_eiz_182 =  empty(len(z))

# CONVERT ABSCISAE FROM EV TO HERTZ:
#----------------------------------------------------------
x_nopeak  = x_nopeak_eV  #*1.5187*10**(15)
x_042_peak= x_042_peak_eV#*1.5187*10**(15)
x_062_peak= x_062_peak_eV#*1.5187*10**(15)
x_082_peak= x_082_peak_eV#*1.5187*10**(15)
x_102_peak= x_102_peak_eV#*1.5187*10**(15)
x_132_peak= x_132_peak_eV#*1.5187*10**(15)
x_182_peak= x_182_peak_eV#*1.5187*10**(15)

## calc difference eps2 and eiz for with and without peak
#-----------------------------------------------------------
#diff_eiz_042 = eiz_042_peak - eiz_nopeak
#diff_eiz_062 = eiz_062_peak - eiz_nopeak
#diff_eiz_082 = eiz_082_peak - eiz_nopeak
#diff_eiz_102 = eiz_102_peak - eiz_nopeak
#diff_eiz_132 = eiz_132_peak - eiz_nopeak
#diff_eiz_182 = eiz_182_peak - eiz_nopeak

listofzeros = np.zeros(len(x_nopeak)) # plot line for y = 0
pl.figure()
for j in range(len(z)):
    eiz_nopeak_arg=empty(len(x_nopeak))
    eiz_042_peak_arg=empty(len(x_042_peak))
    eiz_062_peak_arg=empty(len(x_062_peak))
    eiz_082_peak_arg=empty(len(x_082_peak))
    eiz_102_peak_arg=empty(len(x_102_peak))
    eiz_132_peak_arg=empty(len(x_132_peak))
    eiz_182_peak_arg=empty(len(x_182_peak))
    for i in range(len(x_nopeak)):
        eiz_nopeak_arg[i]=x_nopeak[i]*y_nopeak[i] / (x_nopeak[i]**2 + z[j]**2)
    eiz_nopeak[j] = 1. + (2./pi) * trapz(eiz_nopeak_arg,x_nopeak)

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

    for q in range(len(x_132_peak)):
        eiz_132_peak_arg[q]=x_132_peak[q]*y_132_peak[q]/(x_132_peak[q]**2 + z[j]**2)
    eiz_132_peak[j] = 1. + (2./pi)*trapz(eiz_132_peak_arg,x_132_peak)    

    for s in range(len(x_182_peak)):
        eiz_182_peak_arg[s]=x_182_peak[s]*y_182_peak[s]/(x_182_peak[s]**2 + z[j]**2)
    eiz_182_peak[j] = 1. + (2./pi)*trapz(eiz_182_peak_arg,x_182_peak)    

pl.plot(x_nopeak,eiz_nopeak_arg)
pl.plot(x_182_peak,eiz_182_peak_arg)
pl.close()

savetxt("130723_data/x_nopeak.txt", x_nopeak)
savetxt("130723_data/x_042_peak.txt", x_042_peak)
savetxt("130723_data/x_062_peak.txt", x_062_peak)
savetxt("130723_data/x_082_peak.txt", x_082_peak)
savetxt("130723_data/x_102_peak.txt", x_102_peak)
savetxt("130723_data/x_132_peak.txt", x_132_peak)
savetxt("130723_data/x_182_peak.txt", x_182_peak)
         
savetxt("130723_data/eiz_nopeak.txt", eiz_nopeak)
savetxt("130723_data/eiz_042_peak.txt", eiz_042_peak)
savetxt("130723_data/eiz_062_peak.txt", eiz_062_peak)
savetxt("130723_data/eiz_082_peak.txt", eiz_082_peak)
savetxt("130723_data/eiz_102_peak.txt", eiz_102_peak)
savetxt("130723_data/eiz_132_peak.txt", eiz_132_peak)
savetxt("130723_data/eiz_182_peak.txt", eiz_182_peak)

pl.figure()
pl.plot(  x_nopeak, y_nopeak  , label =r'ab initio' )  
pl.plot(x_042_peak, y_042_peak, label =r'$\omega_{0} = 0.638$')
pl.plot(x_062_peak, y_062_peak, label =r'$\omega_{0} = 0.942$')
pl.plot(x_082_peak, y_082_peak, label =r'$\omega_{0} = 1.245$')
pl.plot(x_102_peak, y_102_peak, label =r'$\omega_{0} = 1.549$')
pl.plot(x_132_peak, y_132_peak, label =r'$\omega_{0} = 2.005$')
pl.plot(x_182_peak, y_182_peak, label =r'$\omega_{0} = 2.764$')
#pl.xlabel(r'$\omega\,\,[sec^{-1}]$', size = 'x-large')
#pl.ylabel(r'$\epsilon"(\omega)$', size = 'x-large')
pl.legend(loc = 'best')
pp.savefig()
pl.show()
pl.close()

pl.figure()
pl.plot(  x_nopeak, y_nopeak  , label =r'ab initio' )  
pl.plot(x_182_peak, y_182_peak, label =r'$\omega_{0} = 2.764$')
#pl.xlabel(r'$\omega\,\,[sec^{-1}]$', size = 'x-large')
#pl.ylabel(r'$\epsilon$"($\omega$)', size = 'x-large')
pl.legend(loc = 'best')
pp.savefig()
pl.show()
pl.close()

pl.figure()
pl.plot(z, eiz_nopeak, label =r'ab initio')  
pl.plot(z, eiz_042_peak, label =r'$\omega_{0} = 0.638$')
pl.plot(z, eiz_062_peak, label =r'$\omega_{0} = 0.942$')
pl.plot(z, eiz_082_peak, label =r'$\omega_{0} = 1.245$')
pl.plot(z, eiz_102_peak, label =r'$\omega_{0} = 1.549$')
pl.plot(z, eiz_132_peak, label =r'$\omega_{0} = 2.005$')
pl.plot(z, eiz_182_peak, label =r'$\omega_{0} = 2.764$')
#pl.xlabel(r'$\xi_{n}\,\,[sec^{-1}]$', size = 'x-large')
#pl.ylabel(r'$\epsilon (i \xi_{n} ) $', size = 'x-large')
pl.legend(loc = 'best')
pp.savefig()
pl.show()
pl.close()

pl.figure()
pl.plot(z, eiz_nopeak, label =r'ab initio')  
pl.plot(z, eiz_182_peak, label =r'$\omega_{0} = 2.764$')
#pl.xlabel(r'$\xi_{n}\,\,[sec^{-1}]$', size = 'x-large')
#pl.ylabel(r'$\epsilon (i \xi_{n} ) $', size = 'x-large')
pl.legend(loc = 'best')
pp.savefig()
pl.show()
pl.close()

diff_eiz_042 = eiz_042_peak - eiz_nopeak
diff_eiz_062 = eiz_062_peak - eiz_nopeak
diff_eiz_082 = eiz_082_peak - eiz_nopeak
diff_eiz_102 = eiz_102_peak - eiz_nopeak
diff_eiz_132 = eiz_132_peak - eiz_nopeak
diff_eiz_182 = eiz_182_peak - eiz_nopeak

pl.figure()
pl.plot(  z, diff_eiz_042, label =r'$\omega_{0} = 0.638$')
pl.plot(  z, diff_eiz_062, label =r'$\omega_{0} = 0.942$')
pl.plot(  z, diff_eiz_082, label =r'$\omega_{0} = 1.245$')
pl.plot(  z, diff_eiz_102, label =r'$\omega_{0} = 1.549$')
pl.plot(  z, diff_eiz_132, label =r'$\omega_{0} = 2.005$')
pl.plot(  z, diff_eiz_182, label =r'$\omega_{0} = 2.764$')
#pl.xlabel(r'$\xi_{n}\,\,[sec^{-1}]$', size = 'x-large')
#pl.ylabel(r'$ \delta\epsilon (i \xi_{n} ) $', size = 'x-large')
pl.legend(loc = 'best')
pl.show()
pp.close()
pl.close()










