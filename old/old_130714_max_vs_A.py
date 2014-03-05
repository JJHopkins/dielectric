#!/usr/bin/python
import matplotlib               
import numpy                    
from pylab import *
from matplotlib import pyplot as pl
from matplotlib import axis as ax
from matplotlib.ticker import MultipleLocator#, FormatStrFormatter
#from scipy import integrate
from matplotlib.backends.backend_pdf import PdfPages
#pp = PdfPages('130712_Hopkins_normed_y_df_DF.pdf')

x_nopeak = numpy.loadtxt('x_nopeak.txt', unpack=True, usecols = [0])
y_nopeak = numpy.loadtxt('a-sio2-eps2.txt', unpack=True, usecols = [1])

x_peak = numpy.loadtxt('x_peak.txt', unpack=True, usecols = [0])
y_peak = numpy.loadtxt('Y24204L.txt', unpack=True, usecols = [1])# L => eps2

x_Re_expt_eV, y_Re_expt = numpy.loadtxt('Y24204K.txt', unpack=True, usecols = [0,1])# K => eps1


eiz_nopeak = numpy.loadtxt('eiz_nopeak.txt', unpack=True, usecols = [0])
eiz_peak = numpy.loadtxt('eiz_peak.txt', unpack=True, usecols = [0])

coeff = (2.41*1e14) # in rad/s
n = arange(0,250)
z = n * coeff
z_y = linspace(0,250,250)
A = 4.0 #[1,4,6]
B = 7.0 #[1,4,6]
C = 10.0 #[1,4,6]

sum_Li2 = numpy.empty(len(z))
dF = numpy.empty(len(z))

sum_Li2_A_dep = numpy.empty(len(z_y))
sum_Li2_yA = numpy.empty(len(z_y))
sum_Li2_yB = numpy.empty(len(z_y))
sum_Li2_yC = numpy.empty(len(z_y))

yA_dep = numpy.empty(len(z_y))
norm_yA_dep = numpy.empty(len(z_y))
yA = numpy.empty(len(z_y))
yB = numpy.empty(len(z_y))
yC = numpy.empty(len(z_y))

A_dep = linspace(1.,11.,300)
#A_dep = [2.,4.,6.,8.,10.,12.]#linspace(5,36,6)
y_0 = numpy.empty(len(A_dep))
sum_dF =  0.0
sum_F =  0.0
sum_F_p =  0.0
diff_sum_F =  0.0
max_yA = numpy.empty(len(A_dep))
max_xA = numpy.empty(len(A_dep))
sum_Li3 = numpy.empty(len(z))
sum_Li3_p = numpy.empty(len(z))
F_3 = numpy.empty(len(z))
F_3_p = numpy.empty(len(z))

colors = ['b','g','r','k','c','m','y']
for r in range(len(A_dep)):
    for t in range(len(z_y)):
        sum_Li2_A_dep[t] = 0.0
        for u in arange(1,101):
            sum_Li2_A_dep[t] += (((A_dep[r]/(A_dep[r]+2.*z_y[t]**2 +2))**2)**u)/(u**2)
        yA_dep[t] = (sum_Li2_A_dep[t]*\
        (A_dep[r]/(A_dep[r]+2.*z_y[t]**2+2))**(-1)*\
        (1-(A_dep[r]/(A_dep[r]+2.*z_y[t]**2+2))**2))\

    #max_xA[r]= z_y[yA_dep.argmax()]
#    pl.plot(z_y,yA_dep/yA_dep[0],  color = colors[r])
#pl.text(8,0.35,r'$A\,=\,2$', color='b')
#pl.text(8,0.30,r'$A\,=\,4$', color='g')
#pl.text(8,0.25,r'$A\,=\,6$', color='r')
#pl.text(8,0.20,r'$A\,=\,8$', color='k')
#pl.text(8,0.15,r'$A\,=\,10$', color='c')
#pl.text(8,0.10,r'$A\,=\,12$', color='m')

#pl.show()
    max_xA[r]= z_y[yA_dep.argmax()]
    #print 2.*10**(16)*max_xA, max_yA  
#pl.xlabel(r'$\frac{\xi}{\omega_{0}}$', size = 'x-large')
#pl.ylabel(r'$ y (\frac{\xi}{\omega_{0}} ) $', size = 'x-large')

pl.figure()
pl.plot(A_dep, max_xA)  
#pl.plot(A_dep,max_yA/yA_dep[0], color = colors[r])
#pl.xlabel(r'$A\,=\,\frac{\2}{\pi}C$', size = 'x-large')
#pl.ylabel(r'$ \frac{\xi}{\omega_{0}} at y_{max}(\frac{xi}{omega_{0}}) $', size = 'x-large')
pl.show()
pl.close()        
#pp.close()
#pl.close()
