#!/usr/bin/python

import numpy as np
from numpy import log,exp
#from scipy import *
import scipy.optimize as opt
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import axis as ax
import matplotlib

import matplotlib               
import numpy                    
from pylab import *
from matplotlib import pyplot as pl
from matplotlib import axis as ax
from matplotlib.ticker import MultipleLocator#, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
# Things will go wrong a lot in the code below
# Ignore most of the warnings
#np.seterr(all='ignore')

Vbar = 1.0
alpha = 0.49
Ns = 9
Nb = 79
#dfs = [0.0, 4.0,9.0]
dfs = np.linspace(0,4,100)#[4.0]
#dfs = [0.0,-5.0, -10.0] #minimized content for poster
#phi_bs = np.linspace(0.00000001, 0.20, 100)
phi_bs = [0.20]#, 0.30]# for poster
phi_water_min = 0.01
phis = np.linspace(0.0000000001, 0.78,100)#-phi_b-phi_water_min, 100)
height= np.zeros(shape=(len(phis),len(dfs)))
# Define eqn that is to be solved numerically
#Variables: p=partition coefficent, phi=vol fract of short, phib=vol frac big
#df=Free energy cost by pore for entering pore, Ns and Nb are number of subunits
# for short and big polymers
def f(p, phi, phib, df):
	""" Implicit definition of P(phi, dF) """
	#return -log(p) - df + (p-1)*phi + \
	#    ( phi*(1-p) + \
	#      5./4*alpha*(phi*p)**(5./4)*(p*phi-9./5) - \
	#      (phib/Nb)-5./4*alpha*(phi+phib)**(5./4)*(phi+phib-9./5) - \
	#      1./2*((1-p*phi)**2 - (1-phi-phib)**2) ) * Ns

	#return -p + exp(- df + (p - 1)*phi + Ns*(phi*(1 - p) + phib + (5./4)*alpha*p**(5./4)*phi**(5./4)*(p*phi - (9./5))- \
	#	phib/Nb - (5./4)*alpha*(phi + phib)**(5./4)*(phi + phib-(9./5))-\
	#	(1./2)*((1 - p*phi)**(2) - (1 - phi - phib)**(2)))) 
	##return -p + exp(- df + Ns*(log((1 - p*phi)/(1 - phi)) + (p - 1)*phi + (9./4)*alpha_tilda*(1 - p**(5./4))*phi**(5./4))
	return - p + exp( - df + Ns*(log((1 - p*phi)/(1 - phi - phib)) + \
		(p - 1)*phi - phib + (9./4)*alpha*((phi + phib)**(5./4) - (p*phi)**(5./4))))

	# Define Osmotic Pressure equation as function of short and big fraction
def osmotic_pressure(s, b):
	#return s/Ns + b/Nb + (5./4)*alpha*(s + b)**(9./4)
	return - log(1 - s - b) + (1 - s - b) - 1 + s/Ns + b/Nb - \
		(1./2)*(1 - (1 - s - b))**2 + (5./4)*alpha* \
		(1 - (1 - s - b))**(9./4)


#Optimization algorithms, if opposite sign use brent, else use newton
def P(phi, phib, df):
	""" Numerically solve for partition coefficient as a
	    function of \phi_s """
	if f(0,phi,phib,df)*f(1,phi,phib,df) < 0:
		#print 'bisect'
		#return opt.brentq(f, 0, 1, args=(phi,df)) # Brent's method
		#return opt.brenth(f, 0, 1, args=(phi,df)) # Brent's method
		return opt.bisect(f, 0, 1, args=(phi,phib,df), maxiter=500) # Bisection method
		#x,r = opt.bisect(f, 0, 1, args=(phi,df), full_output=True) # Bisection method
		#print r.iterations
		#return x
	else:
		#print 'newton'
		return opt.newton(f, 1.0, args=(phi,phib,df), maxiter=5000) # Newton-Raphson
		#return opt.bisect(f, 0, 1, args=(phi,phib,df)) # Bisection method
		# linestyles for different phi_big, colors for different deltaf's

lts = [':','--','-']
#lts = ['--','-']#for poster
markers = ['x', '^', 'd']
#colors = ['k','r','b','g','c','m','y']

for i,phi_b in enumerate(phi_bs):
	#phis = np.linspace(0.0, 1.0-phi_b-phi_water_min, 100)
	for k,phi in enumerate(phis): 
    #	pl.subplot(len(phi_bs),1,i+1)
    #	for df in dfs:
    	    for j,df in enumerate(dfs):
    	    	try: ps = P(phi,phi_b,df)
    	    	#try: ps = [P(phi,phi_b,df) for phi in phis]
    	        except: continue
    	    	#resids = [abs(f(p,phi,phi_b,df)) for p,phi in zip(ps, phis)]
    	    	#max_resid = max(resids)
    	    	#print 'Largest residual for df=%f: %e' % (df, max_resid)
    	    	if j==0:
    	    		labels=r'$\phi_{b} = %.2f$' % phi_b
    	    	else:
    	    		labels=None
    	    	if i==2:
    	    		labelss=r'$\Delta f = %.1f$' % df
    	    	else:
    	    		labelss=None
    	    	#pps = [phi*P(phi,phi_b,df) for phi in phis]
    	    	print 'Last phi tried for df=%f: phis=%.2f: phib=%.1f: total phi=%.2f' % (df, phi,phi_b,phi+phi_b)
    	    	#for k, phi in enumerate(phis):
    	    	#op = [osmotic_pressure(phi,phi_b) for phi in phis]
    	    	op = osmotic_pressure(phi,phi_b)
    	    	#	op[k] = osmotic_pressure(phi,phi_b)
    	    	height[k,j] = ps
    	    	#pl.plot(op,ps,'k', marker= markers[j],linestyle = lts[i], linewidth= 0.5, markersize= 4, markevery= 7, markerfacecolor= 'w')#label=labels)#color = colors[j],
    		#pl.plot(-2.0,-2.0,'k',linestyle = lts[i], linewidth= 1.0,label=labels)#color = colors[j],
    		#pl.plot(-1.0,-1.0,'k',marker= markers[j],linestyle = 'none', markersize= 6,markerfacecolor= 'w',label=labelss)#color = colors[j],
    		#pl.text(0.9,1.95,r'$\Delta F = %.2f$' % j(0), color='r')
    		#pl.text(0.9,1.85,r'$\Delta F = %.2f$' % j(1), color='g')
    		#pl.text(0.9,1.75,r'$\Delta F = %.2f$' % j(2), color='b')
    		#pl.text(0.9,1.65,r'$\Delta F = %.2f$' % j(3), color= 'k')
X,Y = meshgrid(phis,dfs)

figure()
contourf(X,Y,height,1000); clim(0,4)#,cmap = hot())
levels = [0.1,0.5,1.0,1.2,1.3,1.4,1.5,2.0]#,0.80,1.00]
CS = contour(X,Y,height,levels, colors = 'k')
man_loc = [(1,1),(2,2),(3,3),(4,4)]
clabel(CS, inline =1,fmt = '%2.1f', fontsize = 18,color = 'k', manual = man_loc)
#contourf(X,Y,height, 1000); clim(0,4)#;cbar = colorbar();
#CS=contour(X,Y,height, levels)#, colors = 'k')
#clabel(CS,fmt = '%2.1f',colors = 'k',fontsize = 18)
xlabel(r'$\mathrm{\Phi_{PEG\,400}\,\,added\, to\, fixed\,\,\Phi_{PEG\,3500}=\,0.20}$', size = 18)#\,=\,\frac{2}{\pi}C$', size = 20)
ylabel(r'$\Delta f$', size = '20')#\,\,\,Energy cost to enter pore$', size = 20)
#cbar.ax.set_ylabel(r'$P\,=\,\frac{\Phi_{400}^{in}}{\Phi_{400}^{out}}$', size = 20)
#cbar.add_lines(CS)
#axes([0,0.78,0,3])

savefig('plots/130815_Hopkins_PPP_contour_phis_df_P.png', dpi = 300)
#cbar.ax.set_ylabel(r'$\frac{\xi}{\omega_{0}}$', size = 24)
#cbar.add_lines(CS)
show()
close()








#	pl.axis([0, 4.0, -0.10, 3.0])

#font = {'size': 22}
#matplotlib.rc('font',**font)

#pl.text(3.5,1.76,r'$\Delta f = 0.0$', color='k')
#pl.text(3.5,1.64,r'$\Delta f = 5.0$', color='r')
#pl.text(3.5,1.52,r'$\Delta f = 10.0$', color='b')
#pl.text(3.5,1.40,r'$\Delta f = 15.0$', color= 'g')

#pl.legend(loc='lower right')
##pl.title('Partition Coefficient vs Osmotic Pressure, Fixed $\Phi_{PEG3500}$(out)')
#pl.xlabel(r'$\tilde {\Pi} ( \phi_{s} , \phi_{b} )$', size = 'x-large')
#pl.ylabel(r'$p  ( \phi_{s}  , \phi_{b} , \Delta f )$', size = 'x-large')
#pl.axis([0.0, 3.85,-0.01, 1.86])
#
#pl.savefig('VaryShort_BW_PC_vs_OP.eps', dpi=600)
#
#
#pl.show()



