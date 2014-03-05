#!/usr/bin/python
import matplotlib               
from numpy import *                    
from pylab import *
from matplotlib import pyplot as pl
from matplotlib import axis as ax
from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('plots/130801_Hopkins_pre_eV_integral_of_one_over_omega.pdf')

## DEFINE FUNCTIONS FOR CALCULATING e(iz)
#---------------------------------------------------------- 
# Matsubara frequencies: z_n at room temp is (2pikbT/hbar)*n (ie coeff*n)
coeff = 0.159 # in eV#(2.41*1e14) # in rad/s
n = arange(0,250)
z = n * coeff

# INPUT DATA FILES:
#----------------------------------------------------------
#     WAI YIM: ab initio data for a-sio2 (no peak)
x_nopeak_eV, y_nopeak = loadtxt('data/Y23605L.txt', unpack=True, usecols = [0,1])
# DD resent this within set of peak shift data (see x_82_peak_eV):
#x_peak_eV, y_peak = loadtxt('Y24204L.txt', unpack=True, usecols = [0,1])# L => eps2

#     DAN DRYDEN: SYNTH WITH PEAK AT DIFF W_0'S
#----------------------------------------------------------
x_042_peak_eV,y_042_peak=loadtxt('data/Y24901L.txt',unpack=True, usecols = [0,1])#peak at 4.2
x_062_peak_eV,y_062_peak=loadtxt('data/Y24902L.txt',unpack=True, usecols = [0,1])#peak at 6.2
x_082_peak_eV,y_082_peak=loadtxt('data/Y24903L.txt',unpack=True, usecols = [0,1])#peak at 8.2
x_102_peak_eV,y_102_peak=loadtxt('data/Y24904L.txt',unpack=True, usecols= [0,1])#peak at 10.2
x_132_peak_eV,y_132_peak=loadtxt('data/Y24905L.txt',unpack=True, usecols= [0,1])#peak at 13.2
x_182_peak_eV,y_182_peak=loadtxt('data/Y24906L.txt',unpack=True, usecols= [0,1])#peak at 18.2


x_s,y_s=loadtxt('data/Y24204S.PRN',unpack=True, usecols = [0,1])#peak at 4.2
x_042_s,y_042_s=loadtxt('data/Y24901S.PRN',unpack=True, usecols = [0,1])#peak at 4.2
x_062_s,y_062_s=loadtxt('data/Y24902S.PRN',unpack=True, usecols = [0,1])#peak at 6.2
x_082_s,y_082_s=loadtxt('data/Y24903S.PRN',unpack=True, usecols = [0,1])#peak at 8.2
x_102_s,y_102_s=loadtxt('data/Y24904S.PRN',unpack=True, usecols= [0,1])#peak at 10.2
x_132_s,y_132_s=loadtxt('data/Y24905S.PRN',unpack=True, usecols= [0,1])#peak at 13.2
x_182_s,y_182_s=loadtxt('data/Y24906S.PRN',unpack=True, usecols= [0,1])#peak at 18.2

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

perc_042  = 100.*(I_042-I_ab)/I_ab
perc_062  = 100.*(I_062-I_ab)/I_ab
perc_082  = 100.*(I_082-I_ab)/I_ab
perc_102  = 100.*(I_102-I_ab)/I_ab
perc_132  = 100.*(I_132-I_ab)/I_ab
perc_182  = 100.*(I_182-I_ab)/I_ab

print ('Sum rule percent difference for 042 = %s ' % perc_042) 
print ('Sum rule percent difference for 062 = %s ' % perc_062) 
print ('Sum rule percent difference for 082 = %s ' % perc_082) 
print ('Sum rule percent difference for 102 = %s ' % perc_102) 
print ('Sum rule percent difference for 132 = %s ' % perc_132) 
print ('Sum rule percent difference for 182 = %s ' % perc_182) 


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
x_nopeak  = x_nopeak_eV  *1.5187*10**(15)
x_042_peak= x_042_peak_eV*1.5187*10**(15)
x_062_peak= x_062_peak_eV*1.5187*10**(15)
x_082_peak= x_082_peak_eV*1.5187*10**(15)
x_102_peak= x_102_peak_eV*1.5187*10**(15)
x_132_peak= x_132_peak_eV*1.5187*10**(15)
x_182_peak= x_182_peak_eV*1.5187*10**(15)


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
        eiz_nopeak_arg[i]=(x_nopeak[i])**(-1)*y_nopeak[i]
    eiz_nopeak[j] = 1. + (2./pi) * trapz(eiz_nopeak_arg,x_nopeak)

    for k in range(len(x_042_peak)):
        eiz_042_peak_arg[k]=(x_042_peak[k])**(-1)*y_042_peak[k]
    eiz_042_peak[j] = 1. + (2./pi)*trapz(eiz_042_peak_arg,x_042_peak)    

    for m in range(len(x_062_peak)):
        eiz_062_peak_arg[m]=(x_062_peak[m])**(-1)*y_062_peak[m]
    eiz_062_peak[j] = 1. + (2./pi)*trapz(eiz_062_peak_arg,x_062_peak)    

    for n in range(len(x_082_peak)):
        eiz_082_peak_arg[n]=(x_082_peak[n])**(-1)*y_082_peak[n]
    eiz_082_peak[j] = 1. + (2./pi)*trapz(eiz_082_peak_arg,x_082_peak)    

    for p in range(len(x_102_peak)):
        eiz_102_peak_arg[p]=(x_102_peak[p])**(-1)*y_102_peak[p]
    eiz_102_peak[j] = 1. + (2./pi)*trapz(eiz_102_peak_arg,x_102_peak)    

    for q in range(len(x_132_peak)):
        eiz_132_peak_arg[q]=(x_132_peak[q])**(-1)*y_132_peak[q]
    eiz_132_peak[j] = 1. + (2./pi)*trapz(eiz_132_peak_arg,x_132_peak)    

    for s in range(len(x_182_peak)):
        eiz_182_peak_arg[s]=(x_182_peak[s])**(-1)*y_182_peak[s]
    eiz_182_peak[j] = 1. + (2./pi)*trapz(eiz_182_peak_arg,x_182_peak)    

pl.plot(  x_nopeak*10**(-16),  eiz_nopeak_arg*10**(16),  label = r'$\rm{ab\,initio}$')
pl.plot(x_042_peak*10**(-16),eiz_042_peak_arg*10**(16),  label = r'$\omega_{0}= 0.64$')
pl.plot(x_062_peak*10**(-16),eiz_062_peak_arg*10**(16),  label = r'$\omega_{0}= 0.94$')
pl.plot(x_082_peak*10**(-16),eiz_082_peak_arg*10**(16),  label = r'$\omega_{0}= 1.25$')
pl.plot(x_102_peak*10**(-16),eiz_102_peak_arg*10**(16),  label = r'$\omega_{0}= 1.55$')
pl.plot(x_132_peak*10**(-16),eiz_132_peak_arg*10**(16),  label = r'$\omega_{0}= 2.01$')
pl.plot(x_182_peak*10**(-16),eiz_182_peak_arg*10**(16),  label = r'$\omega_{0}= 2.76$')

pl.xlabel(r'$[sec^{-1}]$', size = 'x-large')
pl.ylabel(r'$\rm{arg}\, = \,\omega^{-1}\rm{Im}(\omega)$', size = 'x-large')
pl.title(r'$\rm{\,Integrand\,of}\,\epsilon(i\xi)\,\rm{for\,\xi\,in\,zero\,limit}$')
pl.legend(loc='best')
pp.savefig()
pl.show()
pl.close()

Arg_nopeak=trapz(eiz_nopeak_arg,x_nopeak)
Arg_042   = trapz(eiz_042_peak_arg,x_042_peak)  
Arg_062   = trapz(eiz_062_peak_arg,x_062_peak)  
Arg_082   = trapz(eiz_082_peak_arg,x_082_peak)  
Arg_102   = trapz(eiz_102_peak_arg,x_102_peak)  
Arg_132   = trapz(eiz_132_peak_arg,x_132_peak)  
Arg_182   = trapz(eiz_182_peak_arg,x_182_peak)  

print ('Arg_nopeak = %s ' % Arg_nopeak)
print ('Arg_042 =    %s ' % Arg_042   )
print ('Arg_062 =    %s ' % Arg_062   )
print ('Arg_082 =    %s ' % Arg_082   )
print ('Arg_102 =    %s ' % Arg_102   )
print ('Arg_132 =    %s ' % Arg_132   )
print ('Arg_182 =    %s ' % Arg_182   )





savetxt("data/x_nopeak_0_limit.txt", x_nopeak)
savetxt("data/x_042_peak_0_limit.txt.txt", x_042_peak)
savetxt("data/x_062_peak_0_limit.txt.txt", x_062_peak)
savetxt("data/x_082_peak_0_limit.txt.txt", x_082_peak)
savetxt("data/x_102_peak_0_limit.txt.txt", x_102_peak)
savetxt("data/x_132_peak_0_limit.txt.txt", x_132_peak)
savetxt("data/x_182_peak_0_limit.txt.txt", x_182_peak)
         
savetxt("data/eiz_nopeak_0_limit.txt.txt", eiz_nopeak)
savetxt("data/eiz_042_peak_0_limit.txt.txt", eiz_042_peak)
savetxt("data/eiz_062_peak_0_limit.txt.txt", eiz_062_peak)
savetxt("data/eiz_082_peak_0_limit.txt.txt", eiz_082_peak)
savetxt("data/eiz_102_peak_0_limit.txt.txt", eiz_102_peak)
savetxt("data/eiz_132_peak_0_limit.txt.txt", eiz_132_peak)
savetxt("data/eiz_182_peak_0_limit.txt.txt", eiz_182_peak)

# in eV units
pl.figure()
pl.plot(  x_nopeak_eV, y_nopeak  , label =r'$ab\, initio$' )  
pl.plot(x_042_peak_eV, y_042_peak, label =r'$\omega_{0} =  4.2$')
pl.plot(x_062_peak_eV, y_062_peak, label =r'$\omega_{0} =  6.2$')
pl.plot(x_082_peak_eV, y_082_peak, label =r'$\omega_{0} =  8.2$')
pl.plot(x_102_peak_eV, y_102_peak, label =r'$\omega_{0} = 10.2$')
pl.plot(x_132_peak_eV, y_132_peak, label =r'$\omega_{0} = 13.2$')
pl.plot(x_182_peak_eV, y_182_peak, label =r'$\omega_{0} = 18.2$')
pl.xlabel(r'$\omega\,\,[\rm{eV}]$', size = 'x-large')
pl.ylabel(r'$\epsilon$"($\omega$)', size = 'x-large')
pl.title(r'$\epsilon$"($\omega$) for ab initio and L.O. peaks at $\omega_{0}$')# for ab initio and \omega_{0}=18.2 eV data$')
pl.legend(loc = 'best')
pp.savefig()
pl.show()
pl.close()
# in inverse seconds units
#pl.figure()
#pl.plot(  x_nopeak*1.5187*10**(15), y_nopeak  , label =r'$ab\, initio$' )  
#pl.plot(x_042_peak*1.5187*10**(15), y_042_peak, label =r'$\omega_{0} = 0.64$')
#pl.plot(x_062_peak*1.5187*10**(15), y_062_peak, label =r'$\omega_{0} = 0.94$')
#pl.plot(x_082_peak*1.5187*10**(15), y_082_peak, label =r'$\omega_{0} = 1.25$')
#pl.plot(x_102_peak*1.5187*10**(15), y_102_peak, label =r'$\omega_{0} = 1.55$')
#pl.plot(x_132_peak*1.5187*10**(15), y_132_peak, label =r'$\omega_{0} = 2.01$')
#pl.plot(x_182_peak*1.5187*10**(15), y_182_peak, label =r'$\omega_{0} = 2.76$')
#pl.xlabel(r'$\omega\,\,[sec^{-1}]$', size = 'x-large')
#pl.ylabel(r'$\epsilon$"($\omega$)', size = 'x-large')
#pl.title(r'$\epsilon$"($\omega$) for ab initio and L.O. peaks at $\omega_{0}$')# for ab initio and \omega_{0}=18.2 eV data$')
#pl.legend(loc = 'best')
#pp.savefig()
#pl.show()
#pl.close()

#pl.figure()
#pl.plot(  x_nopeak, y_nopeak  , label =r'ab initio' )  
#pl.plot(x_182_peak, y_182_peak, label =r'$\omega_{0} = 18.2 eV$')
#pl.xlabel(r'$\omega\,\,[sec^{-1}]$', size = 'x-large')
#pl.ylabel(r'$\epsilon$"($\omega$)', size = 'x-large')
#pl.title(r'$\epsilon$"($\omega$) for ab initio and L.O. peak at $\omega_{0}$ = 18.2 eV')# for ab initio and \omega_{0}=18.2 eV data$')
#pl.legend(loc = 'best')
#pp.savefig()
#pl.show()
#pl.close()

#pl.figure()
#pl.plot(z, eiz_nopeak, label =r'$ab\, initio$')  
#pl.plot(z, eiz_042_peak, label =r'$\omega_{0} = 0.64$') 
#pl.plot(z, eiz_062_peak, label =r'$\omega_{0} = 0.94$') 
#pl.plot(z, eiz_082_peak, label =r'$\omega_{0} = 1.25$') 
#pl.plot(z, eiz_102_peak, label =r'$\omega_{0} = 1.55$')
#pl.plot(z, eiz_132_peak, label =r'$\omega_{0} = 2.01$')
#pl.plot(z, eiz_182_peak, label =r'$\omega_{0} = 2.76$')
#pl.xlabel(r'$\xi_{n}\,\,[sec^{-1}]$', size = 'x-large')
#pl.ylabel(r'$\epsilon (i \xi_{n} ) $', size = 'x-large')
#pl.title(r'$\epsilon(i\xi)\, \rm{for\, ab\, initio\, and\, L.O.\, peaks \,at}\, \omega_{0}$')
#pl.legend(loc = 'best')
#pp.savefig()
#pl.show()
#pl.close()

#pl.figure()
#pl.plot(z, eiz_nopeak, label =r'ab initio')  
#pl.plot(z, eiz_182_peak, label =r'$\omega_{0} = 18.2 eV')#2.764$')
#pl.xlabel(r'$\xi_{n}\,\,[sec^{-1}]$', size = 'x-large')
#pl.ylabel(r'$\epsilon (i \xi_{n} ) $', size = 'x-large')
#pl.title(r'$\epsilon(i\xi)$ for ab initio and L.O. peak at $\omega_{0}$ = 18.2 eV$')
#pl.legend(loc = 'best')
#pp.savefig()
#pl.show()
#pl.close()


## calc difference eps2 and eiz for with and without peak
#-----------------------------------------------------------
#diff_eiz_042 = eiz_042_peak - eiz_nopeak
#diff_eiz_062 = eiz_062_peak - eiz_nopeak
#diff_eiz_082 = eiz_082_peak - eiz_nopeak
#diff_eiz_102 = eiz_102_peak - eiz_nopeak
#diff_eiz_132 = eiz_132_peak - eiz_nopeak
#diff_eiz_182 = eiz_182_peak - eiz_nopeak
#
#pl.figure()
#pl.plot(  z/4.2 , diff_eiz_042, color='g', label =r'$\omega_{0} = 0.64$')
#pl.plot(  z/6.2 , diff_eiz_062, color='r', label =r'$\omega_{0} = 0.94$')
#pl.plot(  z/8.2 , diff_eiz_082, color='c', label =r'$\omega_{0} = 1.25$')
#pl.plot(  z/10.2, diff_eiz_102, color='m', label =r'$\omega_{0} = 1.55$')
#pl.plot(  z/13.2, diff_eiz_132, color='y', label =r'$\omega_{0} = 2.01$')
#pl.plot(  z/18.2, diff_eiz_182, color='k', label =r'$\omega_{0} = 2.76$')
#pl.xlabel(r'$\frac{\xi_{n}}{\omega_{0}}$', size = 'x-large')
#pl.ylabel(r'$ \delta\epsilon (i \xi_{n} ) $', size = 'x-large')
#pl.title(r'Difference between $\epsilon^{synth}(i\xi)$ and $\epsilon^{ab\,int}(i\xi)$')
#pl.legend(loc = 'best')
#pp.savefig()
#pl.show()
#pl.close()
#
#pl.figure()
#pl.plot(    x_s*1.5187*10**(15),y_s,  label = r'$ab\,initio$')
#pl.plot(x_042_s*1.5187*10**(15),y_042_s,  label = r'$\omega_0= 0.64$')
#pl.plot(x_062_s*1.5187*10**(15),y_062_s,  label = r'$\omega_0= 0.94$')
#pl.plot(x_082_s*1.5187*10**(15),y_082_s,  label = r'$\omega_0= 1.25$')
#pl.plot(x_102_s*1.5187*10**(15),y_102_s,  label = r'$\omega_0= 1.55$')
#pl.plot(x_132_s*1.5187*10**(15),y_132_s,  label = r'$\omega_0= 2.01$')
#pl.plot(x_182_s*1.5187*10**(15),y_182_s,  label = r'$\omega_0= 2.76$')
#pl.xlabel(r'$\omega\,\,[sec^{-1}]$')#, size = 'x-large')
#pl.ylabel(r'Strength per unit volume (45.3 A$^{3}$)')#, size = 'x-large')
#pl.title('Sum Rule')
##pl.axis([0.0,22.0,0.0,6.0])
#pl.legend(loc='best')
#pp.savefig()
#pl.show()
#pl.close()


#pl.figure()
#pl.plot(    x_s,y_s,  label = r'$ab\,initio$')
#pl.plot(x_042_s/4.2 ,y_042_s,  label = r'$\omega_0 = 0.64$')
#pl.plot(x_062_s/6.2 ,y_062_s,  label = r'$\omega_0 = 0.94$')
#pl.plot(x_082_s/8.2 ,y_082_s,  label = r'$\omega_0 = 1.25$')
#pl.plot(x_102_s/10.2,y_102_s,  label = r'$\omega_0 = 1.55$')
#pl.plot(x_132_s/13.2,y_132_s,  label = r'$\omega_0 = 2.01$')
#pl.plot(x_182_s/18.2,y_182_s,  label = r'$\omega_0 = 2.76$')
#pl.xlabel(r'$\frac{\xi_{n}}{\omega_{0}}$')#, size = 'x-large')
#pl.ylabel(r'Strength per unit volume (45.3 A$^{3}$)')#, size = 'x-large')
#pl.title('Sum Rule')
##pl.axis([0.0,22.0,0.0,6.0])
#pl.legend(loc='best')
#pp.savefig()
#pl.show()
pl.close()
pp.close()










