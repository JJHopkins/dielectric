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

z_y = linspace(0,1.2,500)# z_y is z_n/w_0 ~= e16/e16
A_dep = linspace(0.000001,12,500)

#z_y = linspace(  0.0000000001,1.2 ,200)# z_y is z_n/w_0 ~= e16/e16
#A_dep = linspace(0.0000000000001,7,200)#(0.000001,10,20)
#rho = linspace(  0.0000000000001,5,200)
rho = [1.0,2.0,3.0,4.0]#linspace(  0.0000000000001,5,300)

sum_Li2_A_dep = numpy.empty(len(z_y))
yA_dep = numpy.empty(len(z_y))

#A_dep, rho = meshgrid(A_dep, rho)
labels = np.empty(len(rho))
figure()
height=zeros(shape=(len(A_dep),len(rho)))
for s in range(len(rho)):    
    for r in range(len(A_dep)):
        #print >> sys.stdout, "on r=%d of %d"%(r,len(A_dep))
#    for s in range(len(rho)):    

        for t in range(len(z_y)):
            sum_Li2_A_dep[t] = 0.0
	 
	    U=arange(1,101,1)
	    sum_Li2_A_dep[t] = sum ( ((((A_dep[r]/(1.+(z_y[t]*rho[s])**2))/(2. + (A_dep[r]/(1. + (z_y[t]*rho[s])**2))))**2)**U)/U**2 )
            
            yA_dep[t] = (sum_Li2_A_dep[t])*\
            ((1. - ((A_dep[r]/(1.+(z_y[t]*rho[s])**2)) / (2.+(A_dep[r]/(1.+(z_y[t]*rho[s])**2))))**2)/\
            ((A_dep[r]/(1.+(z_y[t]*rho[s])**2)) / (2.+(A_dep[r]/(1.+(z_y[t]*rho[s])**2)))))*\
            (1./(1. + z_y[t]**2))/(1. + A_dep[r]/(1+(z_y[t]*rho[s])**2)) #  NB this line is different, RP wrote numerator is 1/(1+x^2)...?

        #labels=r'$r\, =\,$',rho[s]#,rho[1],rho[2]


        height[r,s] = z_y[yA_dep==max(yA_dep)] 

X,Y=meshgrid(A_dep,rho)
pl.plot(z_y,height)#, label = labels)
pl.legend(loc = 'upper left')
savetxt('data/131121_fig2_fix0.txt',height)
#
pl.figure()
contourf(X,Y,height,1000,cmap = hot())
CS = contour(X,Y,height, colors = 'k')
#man_loc = [(1,1),(2,2),(3,3),(4,4)]
#clabel(CS, inline =1,fmt = '%2.1f', fontsize = 18,color = 'k', manual = man_loc)
#xlabel(r'$\frac{2}{\pi}C$', size = 21)
#ylabel(r'$r$', size = 24)

pl.xlabel(r'$\varepsilon(0)-1$', size = 'x-large')
pl.ylabel(r'MAX$[\frac{\zeta}{\omega_{0}}] $', size = 'x-large')
#savefig('plots/130815_Hopkins_contour_deltaA_r_xiw0.png', dpi = 300)
show()
close()
