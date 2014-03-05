#!/usr/bin/python
import matplotlib               
import numpy                    
from pylab import *
from matplotlib import pyplot as pl

x = arange(0, 1, 0.1)
y = arange(0, 1, 0.1)
fig = pl.figure()
ax = fig.add_axes([0.1,.1,.8,.8])
ax.plot(x,y)
pl.title('main stuff with long title')


ax_inset = fig.add_axes([0.5,0.5,0.2,0.2])
ax_inset.plot(x,y, label = 'legend')#inset)#, size = 'small')
pl.legend()
pl.title('inset stuff with title', size = 'small')
pl.tick_params(labelsize = 'small')
pl.show()


pl.close()


