#!/usr/bin/python
import matplotlib               
import numpy                    
from pylab import *
from matplotlib import pyplot as pl
from matplotlib import axis as ax
from matplotlib.ticker import MultipleLocator#, FormatStrFormatter
#from scipy import integrate
from matplotlib.backends.backend_pdf import PdfPages
#pp = PdfPages('plots/130814_Hopkins_collagen_eV.pdf')

x_t,y_t= numpy.loadtxt('CollEt.csv', delimiter = ',', unpack=True, usecols = [0,1])
x_x,y_x= numpy.loadtxt('CollEx.csv', delimiter = ',', unpack=True, usecols = [0,1])
x_y,y_y= numpy.loadtxt('CollEy.csv', delimiter = ',', unpack=True, usecols = [0,1])
x_z,y_z= numpy.loadtxt('CollEz.csv', delimiter = ',', unpack=True, usecols = [0,1])


pl.figure()
subplot(111, axisbg ='#F5FFF5')
pl.plot(x_t,y_t, linewidth = '2.0')
pl.plot(x_x,y_x, linewidth = '2.0')
pl.plot(x_y,y_y, linewidth = '2.0')
pl.plot(x_z,y_z, linewidth = '2.0')
pl.xlabel(r'$\hbar\omega\,\,[eV]$',size = 24)#\,=\,\frac{2}{\pi}C
pl.ylabel(r'$\epsilon$"($\omega)$',size = 21)
pl.savefig('130814_Hopkins_hbar_eV_all_collegen_plots_for_DOE.pdf')
#pl.savefig('130813_all_collegen_plots_for_DOE.jpg')


pl.figure()
subplot(111, axisbg ='#F5FFF5')
pl.plot(x_x,y_x,color ='r', linewidth = '2.0')
pl.xlabel(r'$\omega\,\,\,[eV]$',size = 24)#\,=\,\frac{2}{\pi}C
pl.ylabel(r'$\epsilon$"($\omega)$',size = 21)
#pl.savefig('130813_x_collegen_plots_for_DOE.pdf')
pl.savefig('130813_x_collegen_plots_for_DOE.jpg')

pl.figure()
subplot(111, axisbg ='#F5FFF5')
pl.plot(x_y,y_y,color ='r', linewidth = '2.0')
pl.xlabel(r'$\omega\,\,\,[eV]$',size = 24)#\,=\,\frac{2}{\pi}C
pl.ylabel(r'$\epsilon$"($\omega)$',size = 21)
#pl.savefig('130813_y_collegen_plots_for_DOE.pdf')
pl.savefig('130813_y_collegen_plots_for_DOE.jpg')

pl.figure()
subplot(111, axisbg ='#F5FFF5')
pl.plot(x_z,y_z,color ='r', linewidth = '2.0')
pl.xlabel(r'$\omega\,\,\,[eV]$',size = 24)#\,=\,\frac{2}{\pi}C
pl.ylabel(r'$\epsilon$"($\omega)$',size = 21)
#pl.savefig('130813_z_collegen_plots_for_DOE.pdf')
pl.savefig('130813_z_collegen_plots_for_DOE.jpg')
pl.show()

pl.close()        

