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


z_y = linspace(  0.0000000001,1.2 ,300)# z_y is z_n/w_0 ~= e16/e16
A_dep = linspace(0.0000000000001,10,300)#(0.000001,10,20)
rho = [0.9,1.0,1.5,2.0]#linspace(  0.0000000000001,5,300)

sum_Li2_A_dep = numpy.empty(len(z_y))
yA_dep = numpy.empty(len(z_y))

#A_dep, rho = meshgrid(A_dep, rho)
height=zeros(shape=(len(A_dep),len(rho)))
Height=zeros(shape=(len(A_dep),len(rho)))
lts = [':','-.','--','-']
pl.figure()
#for r in range(len(A_dep)):
#    print >> sys.stdout, "on r=%d of %d"%(r,len(A_dep))
for s in range(len(rho)):    

    for r in range(len(A_dep)):
        for t in range(len(z_y)):
            sum_Li2_A_dep[t] = 0.0
	    #A_dep, rho = meshgrid(A_dep, rho)
	 
            ###for u in arange(1,50):	
	    U=arange(1,101,1)
		###sum_Li2_A_dep[t] += ((((A_dep[r]/(1.+(z_y[t]*rho[s])**2))/(2. + (A_dep[r]/(1. + (z_y[t]*rho[s])**2))))**2)**u)/u**2 
	    sum_Li2_A_dep[t] = sum ( ((((A_dep[r]/(1.+(z_y[t]*rho[s])**2))/(2. + (A_dep[r]/(1. + (z_y[t]*rho[s])**2))))**2)**U)/U**2 )
                #sum_Li2_A_dep[t] += (((A_dep[r]/(A_dep[r]+2.*z_y[t]**2 +2))**2)**u)/(u**2)
            #print sum_Li2_A_dep[t]
            #yA_dep[t] = sum_Li2_A_dep[t]*\
            #(A_dep[r]/(A_dep[r]+2.*z_y[p]**2+2))**(-1)*\
            #(1-(A_dep[r]/(A_dep[r]+2.*z_y[t]**2+2))**2)   
            
            yA_dep[t] = (sum_Li2_A_dep[t])*\
            ((1. - ((A_dep[r]/(1.+(z_y[t]*rho[s])**2)) / (2.+(A_dep[r]/(1.+(z_y[t]*rho[s])**2))))**2)/\
            ((A_dep[r]/(1.+(z_y[t]*rho[s])**2)) / (2.+(A_dep[r]/(1.+(z_y[t]*rho[s])**2)))))*\
            (1./(1. + z_y[t]**2))/(1. + A_dep[r]/(1+(z_y[t]*rho[s])**2)) #  NB this line is different, RP wrote numerator is 1/(1+x^2)...?
          
        height[r] = z_y[yA_dep==max(yA_dep)] 
        Height[r,s] = z_y[yA_dep==max(yA_dep)] 
    
    savetxt('data/131128_inset_MRS_fig.txt', Height)
    pl.plot(A_dep, height,color = 'k', linestyle = lts[s])
#    pl.plot(z_y, yA_dep)
X,Y=meshgrid(A_dep,rho)
#figure()
#pl.plot(A_dep, height,color = 'k', linestyle = lts[s])
pl.show()
close()



