#!/usr/bin/python
import matplotlib               
from numpy import *                    
from pylab import *
from matplotlib import pyplot as pl
from matplotlib import axis as ax
#from matplotlib.backends.backend_pdf import PdfPages
#pp = PdfPages('130709_Hopkins_n_from_0_F_dF.pdf')

# wai yim chings ab initio data for a-sio2 (no peak)
x_nopeak_eV, y_nopeak = loadtxt('DATA/Y23605L_np.PRN', unpack=True, usecols = [0,1])
x_peak_eV, y_peak     = loadtxt('DATA/Y24903L_pk.PRN', unpack=True, usecols = [0,1])# L => eps2

x_nopeak = empty(len(x_nopeak_eV))
x_peak   = empty(len(x_peak_eV))

x_nopeak = x_nopeak_eV# * 1.5187*10**(15)
x_peak   =  x_peak_eV #* 1.5187*10**(15)

## DEFINE FUNCTIONS FOR CALCULATING e(iz)
#------------------------------------------------------------- 
# Matsubara frequencies: z_n at room temp is (2pikbT/hbar)*n (ie coeff*n)
coeff = 0.159 # in eV #(2.41*1e14) # in rad/s

n = arange(0,500)
z = n * coeff

eiz_nopeak = empty(len(z))
eiz_peak = empty(len(z))

for j in range(len(z)):
    eiz_nopeak_arg=empty(len(x_nopeak))
    eiz_peak_arg=empty(len(x_peak))

    for i in range(len(x_nopeak)):
        eiz_nopeak_arg[i]=x_nopeak[i]*y_nopeak[i] / (x_nopeak[i]**2 + z[j]**2)
    eiz_nopeak[j] = 1 + (2./pi) * trapz(eiz_nopeak_arg,x_nopeak)

    for k in range(len(x_peak)):
        eiz_peak_arg[k]=x_peak[k]*y_peak[k] / (x_peak[k]**2 + z[j]**2)
	eiz_peak[j] = 1 + (2./pi) * trapz(eiz_peak_arg,x_peak)    

savetxt("DATA/130710_eiz_output_x_nopeak.txt", x_nopeak)
savetxt("DATA/130710_eiz_output_x_peak.txt", x_peak)

savetxt("DATA/130710_eiz_output_eiz_nopeak.txt", eiz_nopeak)
savetxt("DATA/130710_eiz_output_eiz_peak.txt", eiz_peak)
