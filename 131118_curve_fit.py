#!/usr/bin/python
"""
    A simple example using scipy curve_fit to fit data from a file.

    The example provided is a fit of Gaussian or Lorentzian functions
    
    To fit your own data, you need to change:
    (1) def func(x,*p) to return the function you are trying to fit,
    (2) the name of the data file read in by numpy.loadtxt,
    (3) the initial p0 values in the scipy.optimize.curve_fit call.
    http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
"""

import matplotlib               
import numpy as np                    # http://numpy.scipy.org/
import scipy                    # http://scipy.org/
from scipy.optimize import curve_fit #, scipy.special, scipy.stats
from matplotlib import pyplot as pl
import scipy.integrate as grt

def f(x,a,b,c):
    return a*np.exp(-b*x) + c

#def f(x,a,b,c):
#    return a*x/b + c
x_simp = np.linspace(0,5,6)
x = np.linspace(0,500,500)

sum_Li3     = np.loadtxt('data/sum_Li3.txt'    , unpack=True, usecols = [0])
sum_Li3_042 = np.loadtxt('data/sum_Li3_042.txt', unpack=True, usecols = [0])
sum_Li3_062 = np.loadtxt('data/sum_Li3_062.txt', unpack=True, usecols = [0])
sum_Li3_082 = np.loadtxt('data/sum_Li3_082.txt', unpack=True, usecols = [0])
sum_Li3_102 = np.loadtxt('data/sum_Li3_102.txt', unpack=True, usecols = [0])
sum_Li3_132 = np.loadtxt('data/sum_Li3_132.txt', unpack=True, usecols = [0])
sum_Li3_182 = np.loadtxt('data/sum_Li3_182.txt', unpack=True, usecols = [0])

dHs = [2.2199921597 
,1.29587958357 
,0.814230743917 
,0.521828550844 
,0.250977778884 
,0.00301659889764] 
centers = [0.638
,0.942
,1.245
,1.549
,2.005
,2.764]

popt_simp,pcov_simp = curve_fit(f,x_simp,dHs)
fit_simp = popt_simp[0]*np.exp(-popt_simp[1]*x_simp) + popt_simp[2]

popt,pcov = curve_fit(f,x,sum_Li3_182)
print 'popt = ', popt
print 'pcov = ', pcov
fit_182 = popt[0]*np.exp(-popt[1]*x) + popt[2]

coeffs2 = np.polyfit(x_simp, dHs,2)
polynomials2 = np.poly1d(coeffs2)
ys2s = polynomials2(x_simp)

coeffs3 = np.polyfit(x_simp, dHs,3)
polynomials3 = np.poly1d(coeffs3)
ys3s = polynomials3(x_simp)
pl.figure()
pl.plot(centers,dHs)
pl.plot(x_simp,ys2s)
pl.plot(x_simp,ys3s)
pl.plot(x_simp,fit_simp)

coeff2 = np.polyfit(x, sum_Li3_182,2)
polynomial2 = np.poly1d(coeff2)
y2s = polynomial2(x)

coeff3 = np.polyfit(x, sum_Li3_182,3)
polynomial3 = np.poly1d(coeff3)
y3s = polynomial3(x)

coeff4_042 = np.polyfit(x, sum_Li3_042,6)
coeff4_062 = np.polyfit(x, sum_Li3_062,6)
coeff4_082 = np.polyfit(x, sum_Li3_082,6)
coeff4_102 = np.polyfit(x, sum_Li3_102,6)
coeff4_132 = np.polyfit(x, sum_Li3_132,6)
coeff4_182 = np.polyfit(x, sum_Li3_182,6)
polynomial4_042 = np.poly1d(coeff4_042)
polynomial4_062 = np.poly1d(coeff4_062)
polynomial4_082 = np.poly1d(coeff4_082)
polynomial4_102 = np.poly1d(coeff4_102)
polynomial4_132 = np.poly1d(coeff4_132)
polynomial4_182 = np.poly1d(coeff4_182)
y4s_042 = polynomial4_042(x)
y4s_062 = polynomial4_062(x)
y4s_082 = polynomial4_082(x)
y4s_102 = polynomial4_102(x)
y4s_132 = polynomial4_132(x)
y4s_182 = polynomial4_182(x)

pl.figure()
pl.plot(x, sum_Li3_042,color= 'c',linestyle='--')
pl.plot(x, sum_Li3_062,color= 'r',linestyle='--')
pl.plot(x, sum_Li3_082,color= 'g',linestyle='--')
pl.plot(x, sum_Li3_102,color= 'm',linestyle='--')
pl.plot(x, sum_Li3_132,color= 'y',linestyle='--')
pl.plot(x, sum_Li3_182,color= 'k',linestyle='--')
#pl.plot(x, fit_182    , color = 'g' )
#pl.plot(x, y2s        , color = 'r' )
#pl.plot(x, y3s        , color = 'c' )

y4s_042[0] =  1./2 *y4s_042[0]
y4s_062[0] =  1./2 *y4s_062[0]
y4s_082[0] =  1./2 *y4s_082[0]
y4s_102[0] =  1./2 *y4s_102[0]
y4s_132[0] =  1./2 *y4s_132[0]
y4s_182[0] =  1./2 *y4s_182[0]

pl.plot(x, y4s_042, color = 'c')
pl.plot(x, y4s_062, color = 'r')
pl.plot(x, y4s_082, color = 'g')
pl.plot(x, y4s_102, color = 'm')
pl.plot(x, y4s_132, color = 'y')
pl.plot(x, y4s_182, color = 'k')

print 'coeff4_042 = ',coeff4_042[4],coeff4_042[5] 
print 'coeff4_062 = ',coeff4_062[4],coeff4_062[5] 
print 'coeff4_082 = ',coeff4_082[4],coeff4_082[5] 
print 'coeff4_102 = ',coeff4_102[4],coeff4_102[5] 
print 'coeff4_132 = ',coeff4_132[4],coeff4_132[5] 
print 'coeff4_182 = ',coeff4_182[4],coeff4_182[5] 

print 'df4_042 = ', sum(y4s_042)/4.2 - 0.7828535903663599
print 'df4_062 = ', sum(y4s_062)/6.2 - 0.7828535903663599#--0.7873
print 'df4_082 = ', sum(y4s_082)/8.2 - 0.7828535903663599#--0.7873
print 'df4_102 = ', sum(y4s_102)/10.2- 0.7828535903663599#--0.7873
print 'df4_132 = ', sum(y4s_132)/13.2- 0.7828535903663599#--0.7873
print 'df4_182 = ', sum(y4s_182)/18.2- 0.7828535903663599#--0.7873
#pl.figure()
##pl.plot(x, sum_Li3_182, color = 'b')
#pl.plot(x, (fit_182-sum_Li3_182)**2    , color = 'g' )
#pl.plot(x, (y2s -sum_Li3_182   )**2    , color = 'r' )
#pl.plot(x, (y3s -sum_Li3_182   )**2    , color = 'c' )
#pl.plot(x, (y4s -sum_Li3_182   )**2    , color = 'k' )

pl.show()

#savetxt("data/sum_Li3.txt"    , sum_Li3    )
#savetxt("data/sum_Li3_042.txt", sum_Li3_042)
#savetxt("data/sum_Li3_062.txt", sum_Li3_062)
#savetxt("data/sum_Li3_082.txt", sum_Li3_082)
#savetxt("data/sum_Li3_102.txt", sum_Li3_102)
#savetxt("data/sum_Li3_132.txt", sum_Li3_132)
#savetxt("data/sum_Li3_182.txt", sum_Li3_182)
