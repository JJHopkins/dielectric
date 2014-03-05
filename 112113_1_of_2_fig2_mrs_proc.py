#!/usr/bin/python
import matplotlib               
import numpy                    
from pylab import *
from matplotlib import pyplot as pl
from matplotlib import axis as ax
from matplotlib.ticker import MultipleLocator#, FormatStrFormatter
from matplotlib.backends.backend_pdf import PdfPages
#pp = PdfPages('plots/130729_Hopkins_ymax_using_dA.pdf')

z_y = linspace(0,1.2,500)# z_y is z_n/w_0 ~= e16/e16
A_dep = linspace(0.000001,12,500)
sum_Li2_A_dep = numpy.empty(len(z_y))
yA_dep = numpy.empty(len(z_y))
yA_dep_nonnorm = numpy.empty(len(z_y))   

max_yA = numpy.empty(len(A_dep))
max_xA = numpy.empty(len(A_dep))

rho = [0.5,1.0,1.5, 2.0]#linspace(  0.0000000000001,5,300)
#rho = [1.0,2.0,3.0,4.0]#linspace(  0.0000000000001,5,300)
pl.figure()
for s in range(len(rho)):
    for r in range(len(A_dep)):
        for t in range(len(z_y)):
            sum_Li2_A_dep[t] = 0.0
            U=arange(1,101,1)
            #for u in arange(1,101):
            #sum_Li2_A_dep[t] = sum ( ((((A_dep[r]/(1.+(z_y[t])**2))/(2. + (A_dep[r]/(1. + (z_y[t])**2))))**2)**U)/U**2 )
    	    sum_Li2_A_dep[t] = sum ( ((((A_dep[r]/(1.+(z_y[t]*rho[s])**2))/(2. + (A_dep[r]/(1. + (z_y[t]*rho[s])**2))))**2)**U)/U**2 )
                #####sum_Li2_A_dep[t] += ((((A_dep[r]/(1.+z_y[t]**2))/(2. + (A_dep[r]/(1. + z_y[t]**2))))**2)**u)/u**2 
                #sum_Li2_A_dep[t] += (((A_dep[r]/(A_dep[r]+2.*z_y[t]**2 +2))**2)**u)/(u**2)
            
            #yA_dep[t] = sum_Li2_A_dep[t]*\
            #(A_dep[r]/(A_dep[r]+2.*z_y[p]**2+2))**(-1)*\
            #(1-(A_dep[r]/(A_dep[r]+2.*z_y[t]**2+2))**2)   
           
        #yA_dep_nonnorm[t] = (sum_Li2_A_dep[t])*\
        yA_dep[t] = (sum_Li2_A_dep[t])*\
        ((1. - ((A_dep[r]/(1.+(z_y[t]*rho[s])**2)) / (2.+(A_dep[r]/(1.+(z_y[t]*rho[s])**2))))**2)/\
        ((A_dep[r]/(1.+(z_y[t]*rho[s])**2)) / (2.+(A_dep[r]/(1.+(z_y[t]*rho[s])**2)))))*\
        (1./(1. + z_y[t]**2))/(1. + A_dep[r]/(1+(z_y[t]*rho[s])**2)) #  NB this line is different, RP wrote numerator is 1/(1+x^2)...?
  
#        yA_dep[t] = yA_dep_nonnorm[t]/yA_dep_nonnorm[0]  

 
#        yA_dep[t] = (sum_Li2_A_dep[t])*\
#        ((1. - ((A_dep[r]/(1.+z_y[t]**2)) / (2.+(A_dep[r]/(1.+z_y[t]**2))))**2)/\
#        ((A_dep[r]/(1.+z_y[t]**2)) / (2.+(A_dep[r]/(1.+z_y[t]**2)))))*\
#        (1./(1.+z_y[t]**2))/(1. + A_dep[r]/(1+z_y[t]**2))
	
    #max_xA[r]= z_y[yA_dep.argmax()]
    	#print max(yA_dep)
    #max_yA[r] = max(yA_dep)
    ymax = zeros(len(A_dep))
    ymax = max(yA_dep)    
    max_xA[r]= z_y[yA_dep==max(yA_dep)]
pl.plot(A_dep, max_xA)  
    #	pl.plot(A_dep, ymax)  
    #print max_xA
#data = zip(A_dep, max_xA)
#savetxt('data/131125_Hopkins_ymax_output.txt', data)

    #pl.figure()
    #pl.plot(A_dep, max_xA)  
#pl.plot(A_dep,max_yA, color = colors[r])
#pl.plot(A_dep,max_yA/yA_dep[0], color = colors[r])
pl.xlabel(r'$A\,=\,\frac{2}{\pi}C$', size = 'x-large')
pl.ylabel(r'MAX$[y(\frac{\xi}{\omega_{0}})] $', size = 'x-large')
pl.axis([0, 10.2, 0.0, 2.0])
#pp.savefig()
pl.show()
pl.close()
#pp.close()


