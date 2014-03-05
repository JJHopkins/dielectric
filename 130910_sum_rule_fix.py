#!/usr/bin/python
import matplotlib               
from numpy import *                    
#from pylab import *
from matplotlib import pyplot as pl
from matplotlib import axis as ax
#from matplotlib.backends.backend_pdf import PdfPages
#pp = PdfPages('plots/130809_Hopkins_oscillator_strengths.pdf')

## DEFINE FUNCTIONS FOR CALCULATING e(iz)
#---------------------------------------------------------- 
#! Matsubara frequencies: z_n at room temp is (2pikbT/hbar)*n (ie coeff*n)
coeff = 0.159 # in eV#(2.41*1e14) # in rad/s
n = arange(0,500)
z = n * coeff

#! INPUT DATA FILES:
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


#x_npk_s,y_npk_s=loadtxt('data/Y24204S.PRN',unpack=True, usecols = [0,1])#peak at 4.2
x_npk_s,y_npk_ns=loadtxt('data/Y423605S_np.txt',unpack=True, usecols = [0,1])#peak at 4.2
#y_npk_s = 45.3*y_npk_ns
y_npk_s = y_npk_ns
#x_npk_s,y_npk_s=loadtxt('data/ab_int.csv',unpack=True, usecols = [0,1])#peak at 4.2
x_042_s,y_042_s=loadtxt('data/Y24901S.PRN',unpack=True, usecols = [0,1])#peak at 4.2
x_062_s,y_062_s=loadtxt('data/Y24902S.PRN',unpack=True, usecols = [0,1])#peak at 6.2
x_082_s,y_082_s=loadtxt('data/Y24903S.PRN',unpack=True, usecols = [0,1])#peak at 8.2
x_102_s,y_102_s=loadtxt('data/Y24904S.PRN',unpack=True, usecols= [0,1])#peak at 10.2
x_132_s,y_132_s=loadtxt('data/Y24905S.PRN',unpack=True, usecols= [0,1])#peak at 13.2
x_182_s,y_182_s=loadtxt('data/Y24906S.PRN',unpack=True, usecols= [0,1])#peak at 18.2


# CONVERT ABSCISAE FROM EV TO HERTZ:
#----------------------------------------------------------
x_nopeak  = x_nopeak_eV  #*1.5187*10**(15)
x_042_peak= x_042_peak_eV#*1.5187*10**(15)
x_062_peak= x_062_peak_eV#*1.5187*10**(15)
x_082_peak= x_082_peak_eV#*1.5187*10**(15)
x_102_peak= x_102_peak_eV#*1.5187*10**(15)
x_132_peak= x_132_peak_eV#*1.5187*10**(15)
x_182_peak= x_182_peak_eV#*1.5187*10**(15)

listofzeros = zeros(len(x_nopeak)) # plot line for y = 0
ylabels = ['0.0','','1.0','','2.0','','3.0','']
y1slabels = ['0.0','','4.0','','8.0','','12.0','']
yslabels =  ['0.0','','4.0','','8.0','','','']

g = pl.figure()
g, (ax00) = pl.subplots(1,1, sharex = True, sharey = True)
ax00.plot(  x_nopeak, y_nopeak  , color = 'b', label =r'$\rm{calculated}$' )  
ax00.plot(x_042_peak, y_042_peak, color = 'c', label =r'$\omega_{0} =  4.2$')
ax00.plot(x_062_peak, y_062_peak, color = 'r', label =r'$\omega_{0} =  6.2$')
ax00.plot(x_082_peak, y_082_peak, color = 'g', label =r'$\omega_{0} =  8.2$')
ax00.plot(x_102_peak, y_102_peak, color = 'm', label =r'$\omega_{0} = 10.2$')
ax00.plot(x_132_peak, y_132_peak, color = 'y', label =r'$\omega_{0} = 13.2$')
ax00.plot(x_182_peak, y_182_peak, color = 'k', label =r'$\omega_{0} = 18.2$')
ax00.plot(x_182_peak, y_182_peak, color = 'k', linestyle = ':', label =r'$\rm{n_{eff}(\omega)}$')
#ax00.plot(x_182_peak, y_182_peak, color = 'k', linestyle = ':', label =r'$\omega$')
pl.xlabel(r'$\hbar\omega\,\,[eV]$', size = 21)
pl.ylabel(r'$\epsilon$"($\omega$)', size = 21)
#pl.legend(loc = 'best')
#pl.title(r'$\epsilon$"($\omega$) for ab initio and L.O. peaks at $\omega_{0}$')# for ab initio and \omega_{0}=18.2 eV data$')
#pl.savefig('plots/130815_Hopkins_shifted_peaks_eps2.png', dpi = 300)

axs00 = ax00.twinx()
#pl.figure()
axs00.plot(x_npk_s,y_npk_s, color = 'b', linestyle = ':')#,  label = r'$ab\,initio$')
axs00.plot(x_042_s,y_042_s, color = 'c', linestyle = ':')#,  label = r'$\omega_0= 0.64$')
axs00.plot(x_062_s,y_062_s, color = 'r', linestyle = ':')#,  label = r'$\omega_0= 0.94$')
axs00.plot(x_082_s,y_082_s, color = 'g', linestyle = ':')#,  label = r'$\omega_0= 1.25$')
axs00.plot(x_102_s,y_102_s, color = 'm', linestyle = ':')#,  label = r'$\omega_0= 1.55$')
axs00.plot(x_132_s,y_132_s, color = 'y', linestyle = ':')#,  label = r'$\omega_0= 2.01$')
axs00.plot(x_182_s,y_182_s, color = 'k', linestyle = ':')#,  label = r'sum rule')#$\omega_0= 2.76$')
pl.ylabel(r'$n_{\rm{eff}}(\omega)$', size = 21)
#pl.legend(loc = 'best')
pl.axis([0,30,0,5])
g.savefig('Fig4a.pdf')
g.savefig('plots/130911_Hopkins_shifted_eps2_sum_rule.png', dpi = 300)
#pl.savefig('plots/130910_Hopkins_shifted_eps2_sum_rule.png', dpi = 300)
#pp.savefig()

#pp.close()
#pl.close()











