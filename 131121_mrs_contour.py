#!/usr/bin/python
import matplotlib               
import numpy                    
from pylab import *
from matplotlib import pyplot as pl
from matplotlib import axis as ax
from matplotlib.ticker import MultipleLocator#, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
#from matplotlib.backends.backend_pdf import PdfPages
#pp = PdfPages('plots/130814_Hopkins_contour_xiMax_of_r_A.pdf')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

c = loadtxt('output/130807_3D_A_dep.txt', unpack=True, usecols = [0])
C = (numpy.pi/2)*c
r = loadtxt('output/130807_3D_A_rho.txt', unpack=True, usecols = [0])
h = loadtxt('output/130807_3D_height.txt')#, unpack=True, usecols = [0,1])

X,Y=meshgrid(C,r)
h = numpy.nan_to_num(h)

fig = figure()
ax = fig.gca(projection = '3d')
ax.text(-7, 6, 0.7, r'$\zeta/\omega_{0}$', zdir = (-1,1,-3), size = 21)
figure()
#contourf(X,Y,h, 1000, cmap = hot())
surf = ax.plot_surface(X,Y,h, rstride = 20, cstride = 20,alpha = 0.2, cmap = cm.gnuplot, linewidth = 0.5)#gray)#coolwarm)#bone)#hot, linewidth = 0.01, antialiased = True, shade = False)# True)#, cmap = hot()
#surf = ax.plot_wireframe(X,Y,h, rstride = 20, cstride = 20,color = 'k')#True)# cmap = hot, shade = True)#,alpha = 0.9, cmap = cm.hot, linewidth = 0.01, antialiased = True, shade = False)# True)#, cmap = hot()
#colorbar(surf)
#cbar.ax.set_ylabel(r'$\frac{\xi}{\omega_{0}}$', size = 24)
#cset = ax.contour(X,Y,h, zdir = 'z', offset = 0, cmap = cm.jet)
#cset = ax.contour(X,Y,h, zdir = 'x', offset = 5, cmap = cm.jet)
#cset = ax.contourf(X,Y,h, zdir = 'y', offset = 6, cmap = cm.jet)# puts plot of max xi vs discrete r values at r=0 plane
#CS = contour(X,Y,h)#, colors = 'k')
#man_loc = [(1,1),(2,2),(3,3),(4,4)]
#clabel(CS, inline =1,fmt = '%1.1f', fontsize = 18,color = 'k', manual = man_loc)
#ax.grid(on = True)
ax.view_init(elev = 19, azim = -112)
#zlabel(r'$\xi/\omega_{0}$', size = 21)
#ylabel(r'$r$', size = 24)
#xlabel(r'$(\epsilon(0) -1)$', size = 24)
#text = Axes.text(self, x, y, s, **kwargs)
#art3d.text_2d_to_3d(text, z, zdir)
#return text

#pl.text(6,0, 0, r'$\xi/\omega_{0}$',size = 21 ,rotation = 'horizontal')
#ax.text(r'$\xi/\omega_{0}$',6,0, 0, size = 21 ,rotation = 'horizontal')
#ax.set_zlabel(r'$\xi/\omega_{0}$',size = 21 ,rotation = 'horizontal' )
ax.set_xlabel(r'$\epsilon(0)-1$', size = 21)
ax.set_ylabel(r'$r$', size = 22)

#ax.set_zlabel(r'$\xi/\omega_{0}$', size = 22)
#ax.axhline(y = 1, color = 'k', linewidth = 2)
#show()
#savefig('Fig1b.pdf')
####savefig('plots/fig1b.eps')
#cbar.ax.set_ylabel(r'$\frac{\xi}{\omega_{0}}$', size = 24)
#cbar.add_lines(CS)
#
pl.figure()
pl.plot()


