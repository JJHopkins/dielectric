#!/usr/bin/python
import matplotlib               
from numpy import *                    
from pylab import *
from matplotlib import pyplot as pl
from matplotlib import axis as ax
from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('plots/130729_Hopkins_my_high_w0_LO_added_output.pdf')



# A lorentzian peak with:
#   Peak height above background : p[0]
#   Central value                : p[1]
#   Full Width at Half Maximum   : p[2]
#return (p[0])/(1.0+((x-p[1])/p[2])**2)

x_data, y_data = loadtxt('data/a-sio2-eps2.txt', unpack=True, usecols = [0,1])
y_233_plus = empty(len(x_data))
y_263_plus = empty(len(x_data))
y_283_plus = empty(len(x_data))
y_316_plus = empty(len(x_data))
y_363_plus = empty(len(x_data))
y_430_plus = empty(len(x_data))

y_233=2.0 /(1.+((x_data-23.3)/0.5)**2)
y_263=2.0 /(1.+((x_data-26.3)/0.5)**2)
y_283=2.0 /(1.+((x_data-28.3)/0.5)**2)
y_316=2.0 /(1.+((x_data-31.6)/0.5)**2)
y_363=2.0 /(1.+((x_data-36.3)/0.5)**2)
y_430=2.0 /(1.+((x_data-43.0)/0.5)**2)

y_233_plus = y_data + y_233 
y_263_plus = y_data + y_263
y_283_plus = y_data + y_283
y_316_plus = y_data + y_316
y_363_plus = y_data + y_363
y_430_plus = y_data + y_430

pl.figure()
pl.plot(x_data,y_data)
pl.plot(x_data,y_233_plus)
pl.plot(x_data,y_263_plus)
pl.plot(x_data,y_283_plus)
pl.plot(x_data,y_316_plus)
pl.plot(x_data,y_363_plus)
pl.plot(x_data,y_430_plus)
pp.savefig()
pl.show()
#y_data[x_data>=xswitch]=y_43[x_data>=xswitch]
#y_data[x_data>=40.93]=y_43[x_data>=40.93]
#pl.plot(x_data,y_data)
#pl.show()
#high_w0_43_data = zip(x_data,y_data)
savetxt("output/my_w0_233.txt", y_233_plus)
savetxt("output/my_w0_263.txt", y_263_plus)
savetxt("output/my_w0_283.txt", y_283_plus)
savetxt("output/my_w0_316.txt", y_316_plus)
savetxt("output/my_w0_363.txt", y_363_plus)
savetxt("output/my_w0_430.txt", y_430_plus)

pl.close()
