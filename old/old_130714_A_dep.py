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

x_nopeak = numpy.loadtxt('data/x_nopeak.txt', unpack=True, usecols = [0])
y_nopeak = numpy.loadtxt('data/a-sio2-eps2.txt', unpack=True, usecols = [1])

x_peak = numpy.loadtxt('data/x_peak.txt', unpack=True, usecols = [0])
y_peak = numpy.loadtxt('data/Y24204L.txt', unpack=True, usecols = [1])# L => eps2

x_Re_expt_eV, y_Re_expt = numpy.loadtxt('data/Y24204L.txt', unpack=True, usecols = [0,1])# not correct file hereK => eps1


eiz_nopeak = numpy.loadtxt('data/eiz_nopeak.txt', unpack=True, usecols = [0])
eiz_peak = numpy.loadtxt('data/eiz_peak.txt', unpack=True, usecols = [0])

coeff = (2.41*1e14) # in rad/s
n = arange(0,250)
z = n * coeff
z_y = linspace(0,10,250)
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

A_dep = linspace(5,36,50)
y_0 = numpy.empty(len(A_dep))
sum_dF =  0.0
sum_F =  0.0
sum_F_p =  0.0
diff_sum_F =  0.0
max_yA = numpy.empty(len(A_dep))
max_yA_non = numpy.empty(len(A_dep))
sum_Li3 = numpy.empty(len(z))
sum_Li3_p = numpy.empty(len(z))
F_3 = numpy.empty(len(z))
F_3_p = numpy.empty(len(z))

for p in range(len(z_y)):
    sum_Li2_yA[p] = 0.0
    sum_Li2_yB[p] = 0.0
    sum_Li2_yC[p] = 0.0
    for q in arange(1,101):
        #sum_Li2_yA[p] += (((A/(1+z_[p]))**q)/(q**2)# +  (((a/(a+2.*z_y[p]**2 +2)))**q)/(q**2))
        sum_Li2_yA[p] += (((A/(A+2.*z_y[p]**2 +2))**2)**q)/(q**2)
        sum_Li2_yB[p] += (((B/(B+2.*z_y[p]**2 +2))**2)**q)/(q**2)
        sum_Li2_yC[p] += (((C/(C+2.*z_y[p]**2 +2))**2)**q)/(q**2)
    yA[p] = sum_Li2_yA[p]*((A/(A+2.*z_y[p]**2 +2))**(-1)*(1 - (A/(A+2.*z_y[p]**2 +2))**2))   
    yB[p] = sum_Li2_yB[p]*((B/(B+2.*z_y[p]**2 +2))**(-1)*(1 - (B/(B+2.*z_y[p]**2 +2))**2))  
    yC[p] = sum_Li2_yC[p]*((C/(C+2.*z_y[p]**2 +2))**(-1)*(1 - (C/(C+2.*z_y[p]**2 +2))**2)) 
pl.figure()
for r in range(len(A_dep)):
    for t in range(len(z_y)):
        sum_Li2_A_dep[t] = 0.0
        for u in arange(1,101):
            sum_Li2_A_dep[t] += (((A_dep[r]/(A_dep[r]+2.*z_y[t]**2 +2))**2)**u)/(u**2)
        
        yA_dep[t] = sum_Li2_A_dep[t]*\
        (A_dep[r]/(A_dep[r]+2.*z_y[p]**2+2))**(-1)*\
        (1-(A_dep[r]/(A_dep[r]+2.*z_y[t]**2+2))**2)   
        y_0[r] = yA_dep[0]
        norm_yA_dep[r] = yA_dep[t]/y_0[r]#A_dep[0]
        max_yA[r] = max(norm_yA_dep)
	#y_norm[r] = yA_dep[r]/y_0[r]
    #print y_0
    #max_yA_non[r] = max(yA_dep)
    #max_yA = max(yA_dep/y_0)
    pl.plot(A_dep, max_yA)
    #pl.plot(z_y,yA_dep/yA_dep[0])
    pl.show
      #  norm_yA_dep[t] = yA_dep[t]/yA_dep[0]
      #  max_yA[t] = max(norm_yA_dep)
#print A_dep,max_yA
    #max_yA[r] = max(yA_dep)
   # print max_yA

#pl.figure()
#pl.plot(A_dep, max_yA)
#pl.show()

#pp.close()
pl.close()