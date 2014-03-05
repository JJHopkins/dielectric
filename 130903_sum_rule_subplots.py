#!/usr/bin/python
import matplotlib               
from numpy import *                    
#from pylab import *
from matplotlib import pyplot as pl
from matplotlib import axis as ax
from matplotlib.backends.backend_pdf import PdfPages
#pp = PdfPages('plots/130809_Hopkins_oscillator_strengths.pdf')

pl.close('all')
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
x_042_peak_eV,y_042_peak=loadtxt('data/Y24901L.txt',unpack=True, usecols = [0,1])#peak at 4.2
x_062_peak_eV,y_062_peak=loadtxt('data/Y24902L.txt',unpack=True, usecols = [0,1])#peak at 6.2
x_082_peak_eV,y_082_peak=loadtxt('data/Y24903L.txt',unpack=True, usecols = [0,1])#peak at 8.2
x_102_peak_eV,y_102_peak=loadtxt('data/Y24904L.txt',unpack=True, usecols= [0,1])#peak at 10.2
x_132_peak_eV,y_132_peak=loadtxt('data/Y24905L.txt',unpack=True, usecols= [0,1])#peak at 13.2
x_182_peak_eV,y_182_peak=loadtxt('data/Y24906L.txt',unpack=True, usecols= [0,1])#peak at 18.2


#x_npk_s,y_npk_s=loadtxt('data/Y24204S.PRN',unpack=True, usecols = [0,1])#peak at 4.2
x_npk_s,y_npk_ns=loadtxt('data/Y423605S_np.txt',unpack=True, usecols = [0,1])#peak at 4.2
y_npk_s = 45.3*y_npk_ns
#x_npk_s,y_npk_s=loadtxt('data/ab_int.csv',unpack=True, usecols = [0,1])#peak at 4.2
x_042_s,y_042_s=loadtxt('data/Y24901S.PRN',unpack=True, usecols = [0,1])#peak at 4.2
x_062_s,y_062_s=loadtxt('data/Y24902S.PRN',unpack=True, usecols = [0,1])#peak at 6.2
x_082_s,y_082_s=loadtxt('data/Y24903S.PRN',unpack=True, usecols = [0,1])#peak at 8.2
x_102_s,y_102_s=loadtxt('data/Y24904S.PRN',unpack=True, usecols= [0,1])#peak at 10.2
x_132_s,y_132_s=loadtxt('data/Y24905S.PRN',unpack=True, usecols= [0,1])#peak at 13.2
x_182_s,y_182_s=loadtxt('data/Y24906S.PRN',unpack=True, usecols= [0,1])#peak at 18.2

I_ab = trapz(y_npk_s,x_npk_s)
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

f = pl.figure()
f, (ax1,ax2,ax3,ax4,ax5,ax6,ax7) = pl.subplots(7,1, sharex = True, sharey = True)
ax1.plot(  x_nopeak, y_nopeak  , color = 'b', linewidth = 1.5)
ax1.set_xlim(0,43)
ax1.set_ylim(0,4)
ax1.set_yticklabels(ylabels)
ax2.plot(x_042_peak, y_042_peak, color = 'c', linewidth = 1.5)
ax2.set_xlim(0,43)
ax2.set_ylim(0,4)
ax3.plot(x_062_peak, y_062_peak, color = 'r', linewidth = 1.5)
ax3.set_xlim(0,43)
ax3.set_ylim(0,4)
ax4.plot(x_082_peak, y_082_peak, color = 'g', linewidth = 1.5)
ax4.set_xlim(0,43)
ax4.set_ylim(0,4)
ax4.set_ylabel(r'$\epsilon$"($\omega$)', size = 24)
ax5.plot(x_102_peak, y_102_peak, color = 'm', linewidth = 1.5)
ax5.set_xlim(0,43)
ax5.set_ylim(0,4)
ax6.plot(x_132_peak, y_132_peak, color = 'y', linewidth = 1.5)
ax6.set_xlim(0,43)
ax6.set_ylim(0,4)
ax7.plot(x_182_peak, y_182_peak, color = 'k', linewidth = 1.5)
ax7.set_xlim(0,43)
ax7.set_ylim(0,4)
ax7.set_xlabel(r'$\hbar\omega\,\,\,[eV]$', size = 24)
#pl.axis([0,6,0,25])

axs1 = ax1.twinx()
axs1.plot(x_npk_s,y_npk_s, color = 'b', linewidth=1.5, linestyle = '--')
axs1.set_xlim(0,43)
axs1.set_yticklabels(y1slabels)
axs2 = ax2.twinx()
axs2.plot(x_042_s,y_042_s, color = 'c', linewidth=1.5, linestyle = '--')
axs2.set_xlim(0,43)
axs2.set_yticklabels(yslabels)
axs3 = ax3.twinx()
axs3.plot(x_062_s,y_062_s, color = 'r', linewidth=1.5, linestyle = '--')
axs3.set_xlim(0,43)
axs3.set_yticklabels(yslabels)
axs4 = ax4.twinx()
axs4.plot(x_082_s,y_082_s, color = 'g', linewidth=1.5, linestyle = '--')
axs4.set_xlim(0,43)
axs4.set_yticklabels(yslabels)
axs4.set_ylabel(r'Strength', size = 21)
axs5 = ax5.twinx()
axs5.plot(x_102_s,y_102_s, color = 'm', linewidth=1.5, linestyle = '--')
axs5.set_xlim(0,43)
axs5.set_yticklabels(yslabels)
axs6 = ax6.twinx()
axs6.plot(x_132_s,y_132_s, color = 'y', linewidth=1.5, linestyle = '--')
axs6.set_xlim(0,43)
axs6.set_yticklabels(yslabels)
axs7 = ax7.twinx()
axs7.plot(x_182_s,y_182_s, color = 'k', linewidth=1.5, linestyle = '--')
axs7.set_xlim(0,43)
axs7.set_yticklabels(yslabels)
f.subplots_adjust(hspace = 0)
#ax.set_xlabel(r'$\hbar\omega\,\,\,[eV]$')
#axs = fig.add_subplot(711)
#axs.set_ylabel('Strength')
f.savefig('plots/130909_shifted_peaks_table_of_sum_rule_and_eps2.pdf')
##pp.savefig()


g = pl.figure()
g, (ax00) = pl.subplots(1,1, sharex = True, sharey = True)
#ax00.plot(  x_nopeak, y_nopeak  , color = 'b', label =r'$\mathrm{calculated}$' )  
ax00.plot(  x_nopeak, y_nopeak  , color = 'b', label ='calculated' )  
ax00.plot(x_042_peak, y_042_peak, color = 'c', label =r'$\omega_{0} =  4.2$')
ax00.plot(x_062_peak, y_062_peak, color = 'r', label =r'$\omega_{0} =  6.2$')
ax00.plot(x_082_peak, y_082_peak, color = 'g', label =r'$\omega_{0} =  8.2$')
ax00.plot(x_102_peak, y_102_peak, color = 'm', label =r'$\omega_{0} = 10.2$')
ax00.plot(x_132_peak, y_132_peak, color = 'y', label =r'$\omega_{0} = 13.2$')
ax00.plot(x_182_peak, y_182_peak, color = 'k', label =r'$\omega_{0} = 18.2$')
#ax00.plot(x_182_peak, y_182_peak, color = 'k', linestyle = ':', label =r'$\rm{n_{eff}(\omega)}$')
ax00.plot(x_182_peak, y_182_peak, color = 'k', linestyle = ':', label =r'$\omega$')
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
#pl.ylabel(r'$n_{eff}(\omega)$', size = 21)
pl.ylabel(r'$\omega$', size = 21)
#pl.legend(loc = 'best')
pl.axis([0,30,0,5])
g.savefig('plots/130910_Hopkins_shifted_eps2_sum_rule.png', dpi = 300)
#pl.savefig('plots/130910_Hopkins_shifted_eps2_sum_rule.png', dpi = 300)
##pp.savefig()

#
pl.figure()
pl.plot(x_npk_s,y_npk_s, color = 'b', label =r'$calculated$' )        #linestyle = ':')#,  label = r'$ab\,initio$')
pl.plot(x_042_s,y_042_s, color = 'c', label =r'$\omega_{0} =  4.2$')  #linestyle = ':')#,  label = r'$\omega_0= 0.64$')
pl.plot(x_062_s,y_062_s, color = 'r', label =r'$\omega_{0} =  6.2$')  #linestyle = ':')#,  label = r'$\omega_0= 0.94$')
pl.plot(x_082_s,y_082_s, color = 'g', label =r'$\omega_{0} =  8.2$')  #linestyle = ':')#,  label = r'$\omega_0= 1.25$')
pl.plot(x_102_s,y_102_s, color = 'm', label =r'$\omega_{0} = 10.2$')  #linestyle = ':')#,  label = r'$\omega_0= 1.55$')
pl.plot(x_132_s,y_132_s, color = 'y', label =r'$\omega_{0} = 13.2$')  #linestyle = ':')#,  label = r'$\omega_0= 2.01$')
pl.plot(x_182_s,y_182_s, color = 'k', label =r'$\omega_{0} = 18.2$')  #linestyle = ':',    label = r'$\omega_0= 2.76$')
pl.legend(loc = 'best')
#pl.axis([0,30,0,5])
pl.xlabel(r'$\hbar\omega\,\,[eV]$', size = 21)
pl.ylabel(r'$Strength$', size = 21)
pl.savefig('plots/130908_Hopkins_sum_rules.pdf')
#pp.savefig()
#
#####################################################################
#Diffferences:

pl.figure()
#pl.plot(  x_nopeak, y_nopeak  , label =r'$ab\, initio$' )  
pl.plot(x_042_peak, y_042_peak-0.977*y_nopeak, color = 'c')#, label =r'$\omega_{0} = 0.64$')
pl.plot(x_062_peak, y_062_peak-0.977*y_nopeak, color = 'r')#, label =r'$\omega_{0} = 0.94$')
pl.plot(x_082_peak, y_082_peak-0.977*y_nopeak, color = 'g')#, label =r'$\omega_{0} = 1.25$')
pl.plot(x_102_peak, y_102_peak-0.977*y_nopeak, color = 'm')#, label =r'$\omega_{0} = 1.55$')
pl.plot(x_132_peak, y_132_peak-0.977*y_nopeak, color = 'y')#, label =r'$\omega_{0} = 2.01$')
pl.plot(x_182_peak, y_182_peak-0.977*y_nopeak, color = 'k', label =r'spectrum - calculated spectrum')# label =r'$\omega_{0} = 2.76$')
#pl.plot(x_npk_s,y_npk_s,  label = r'$ab\,initio$')
pl.plot(x_042_s,y_042_s-y_npk_s, color = 'c', linestyle = ':')#,  label = r'$\omega_0= 0.64$')
pl.plot(x_062_s,y_062_s-y_npk_s, color = 'r', linestyle = ':')#,  label = r'$\omega_0= 0.94$')
pl.plot(x_082_s,y_082_s-y_npk_s, color = 'g', linestyle = ':')#,  label = r'$\omega_0= 1.25$')
pl.plot(x_102_s,y_102_s-y_npk_s, color = 'm', linestyle = ':')#,  label = r'$\omega_0= 1.55$')
pl.plot(x_132_s,y_132_s-y_npk_s, color = 'y', linestyle = ':')#,  label = r'$\omega_0= 2.01$')
pl.plot(x_182_s,y_182_s-y_npk_s, color = 'k', linestyle = ':', label =r'strength - strength of calculated')#  label = r'$\omega_0= 2.76$')
#pl.axis([0,23,0,6])
pl.legend(loc = 'best')
pl.xlabel(r'$\hbar\omega\,\,[eV]$', size = 21)
#pl.ylabel(r'$\epsilon$"($\omega$)', size = 21)
pl.axis([0,30,-0.5,4.0])
pl.savefig('plots/130905_shifted_peaks_differences_for_sum_rule_and_eps2.pdf')
#pp.savefig()

pl.figure()
pl.plot(x_082_peak,            y_nopeak, color = 'b', label =r'$calc\, spectrum$')
pl.plot(x_082_peak, y_082_peak-0.977*y_nopeak, color = 'g', label =r'$added\, peak$')
pl.plot(x_082_s,                y_082_s, color = 'g', linestyle = ':', marker = 'x', markevery = 100, label = r'$strength\, of\, synth$')
pl.plot(x_082_s,                y_npk_s, color = 'b', linestyle = ':',  label = r'$strength\, of\, calc$')
pl.plot(x_082_s,        y_082_s-y_npk_s, color = 'g', linestyle = '--',  label = r'$strength\, difference$')
pl.legend(loc = 'upper left')
pl.xlabel(r'$\hbar\omega\,\,[eV]$', size = 21)
#pl.ylabel(r'$\epsilon$"($\omega$)', size = 21)
pl.axis([0,30,-0.5,4.0])
pl.savefig('plots/130906_ab_initio_synth_sum_rule question.pdf')

#pp.savefig()
pl.figure()
pl.plot(x_042_s,y_042_s-y_npk_s, color = 'c', label = r'$\omega_0=  4.2$')
pl.plot(x_062_s,y_062_s-y_npk_s, color = 'r', label = r'$\omega_0=  6.2$')
pl.plot(x_082_s,y_082_s-y_npk_s, color = 'g', label = r'$\omega_0=  8.2$')
pl.plot(x_102_s,y_102_s-y_npk_s, color = 'm', label = r'$\omega_0= 10.2$')
pl.plot(x_132_s,y_132_s-y_npk_s, color = 'y', label = r'$\omega_0= 13.2$')
pl.plot(x_182_s,y_182_s-y_npk_s, color = 'k', label = r'$\omega_0= 18.2$')
#pl.axis([0,23,0,6])
pl.legend(loc = 'best')
pl.xlabel(r'$\hbar\omega\,\,[eV]$', size = 21)
pl.ylabel(r'$\Delta Strength$', size = 21)
#pl.savefig('plots/130906_shifted_peaks_sume_rule_differences.pdf')
pl.savefig('plots/130911_Hopkins_shifted_eps2_and_sum_rule.png', dpi = 300)
#pp.savefig()
##pl.figure()
#pl.plot(x_npk_s*1.5187*10**(15),y_npk_s,  label = r'$ab\,initio$')
#pl.plot(x_042_s*1.5187*10**(15),y_042_s,  label = r'$\omega_0= 0.64$')
#pl.plot(x_062_s*1.5187*10**(15),y_062_s,  label = r'$\omega_0= 0.94$')
#pl.plot(x_082_s*1.5187*10**(15),y_082_s,  label = r'$\omega_0= 1.25$')
#pl.plot(x_102_s*1.5187*10**(15),y_102_s,  label = r'$\omega_0= 1.55$')
#pl.plot(x_132_s*1.5187*10**(15),y_132_s,  label = r'$\omega_0= 2.01$')
#pl.plot(x_182_s*1.5187*10**(15),y_182_s,  label = r'$\omega_0= 2.76$')

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
##pp.savefig()
#pl.show()
#pl.show()
#pp.close()
pl.close()










