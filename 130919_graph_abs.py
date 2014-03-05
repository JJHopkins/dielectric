#!/usr/bin/python
import matplotlib               
import numpy                    
from pylab import *
from matplotlib import axis as ax
from matplotlib.ticker import MultipleLocator#, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
#from matplotlib.backends.backend_pdf import PdfPages
#pp = PdfPages('plots/130814_Hopkins_contour_xiMax_of_r_A.pdf')

x_nopeak_eV  ,y_nopeak = loadtxt('DATA/Y23605L_np.PRN', unpack=True, usecols = [0,1])
x_042_peak_eV,y_042_peak=loadtxt('data/Y24901L.txt',unpack=True, usecols = [0,1])#peak at 4.2
x_062_peak_eV,y_062_peak=loadtxt('data/Y24902L.txt',unpack=True, usecols = [0,1])#peak at 6.2
x_082_peak_eV,y_082_peak=loadtxt('data/Y24903L.txt',unpack=True, usecols = [0,1])#peak at 8.2
x_102_peak_eV,y_102_peak=loadtxt('data/Y24904L.txt',unpack=True, usecols= [0,1])#peak at 10.2
x_132_peak_eV,y_132_peak=loadtxt('data/Y24905L.txt',unpack=True, usecols= [0,1])#peak at 13.2
x_182_peak_eV,y_182_peak=loadtxt('data/Y24906L.txt',unpack=True, usecols= [0,1])#peak at 18.2

#om_000 = zeros(len(x_nopeak_eV))
#om_npk = om_000 + 1.0
#c = loadtxt('output/130807_3D_A_dep.txt', unpack=True, usecols = [0])
#C = (numpy.pi/2)*c
#r = loadtxt('output/130807_3D_A_rho.txt', unpack=True, usecols = [0])
#h = loadtxt('output/130807_3D_height.txt')#, unpack=True, usecols = [0,1])

X,Y=meshgrid(x_nopeak_eV, y_nopeak)
h_000 = zeros(len(x_nopeak_eV), len(x_nopeak_eV)) 
h_111 = zeros(len(x_nopeak_eV)) + 1.0

fig = figure()
ax = fig.gca(projection = '3d')
#ax.text(-7, 6, 0.7, r'$\zeta/\omega_{0}$', zdir = (-1,1,-3), size = 21)
surf = ax.plot_surface(X,Y,h_npk)#, rstride = 1, cstride = 1,alpha = 0.5, cmap = cm.jet, linewidth = 0.01, antialiased = True, shade = False)# True)#, cmap = hot()
#cset = ax.contour(X,Y,h, zdir = 'z', offset = 0, cmap = cm.jet)
#cset = ax.contour(X,Y,h, zdir = 'y', offset = 0, cmap = cm.binary)# puts plot of max xi vs discrete r values at r=0 plane
#ax.view_init(elev = 34, azim = -135)

#ax.set_xlabel(r'$\epsilon(0)-1$', size = 21)
#ax.set_ylabel(r'$r$', size = 22)

savefig('130919_abs.pdf')

#savefig('plots/130911_Hopkins_contour_deltaA_r_xiw0.png', dpi = 300)



