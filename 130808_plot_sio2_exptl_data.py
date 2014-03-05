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
x,y = loadtxt('data/asio2_exptl.txt', unpack=True, usecols = [0,1])


pl.figure()
#subplot(111, axisbg = 'k')# '#CCCFCC')
plot(x,y,color = 'r')#, linewidth='3.5',linestyle = '-')
pl.xlabel(r'$\hbar\omega\,\,[eV]$', size = 21)
pl.ylabel(r'$\epsilon$"($\omega$)', size = 21)
pl.axis([0,45,-0.2,4.5])
#pl.savefig('plots/130814_aSiO2_experimental_spectrum.pdf')
pl.savefig('plots/130815_Hopkins_aSiO2_experimental_spectrum.png', dpi = 300)
show()
pl.close()








