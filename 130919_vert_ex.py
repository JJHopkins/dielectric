from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
from matplotlib.colors import colorConverter
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

fig = plt.figure()
ax = fig.gca(projection='3d')
cc = lambda arg: colorConverter.to_rgba(arg, alpha=0.9)

x            ,y        = np.loadtxt('data/asio2_exptl.txt',unpack=True, usecols = [0,1])
x_nopeak_eV  ,y_nopeak = np.loadtxt('DATA/Y23605L_np.PRN', unpack=True, usecols = [0,1])
x_042_peak_eV,y_042_peak=np.loadtxt('data/Y24901L.txt', unpack=True, usecols = [0,1])#peak at 4.2
x_062_peak_eV,y_062_peak=np.loadtxt('data/Y24902L.txt', unpack=True, usecols = [0,1])#peak at 6.2
x_082_peak_eV,y_082_peak=np.loadtxt('data/Y24903L.txt', unpack=True, usecols = [0,1])#peak at 8.2
x_102_peak_eV,y_102_peak=np.loadtxt('data/Y24904L.txt', unpack=True, usecols= [0,1])#peak at 10.2
x_132_peak_eV,y_132_peak=np.loadtxt('data/Y24905L.txt', unpack=True, usecols= [0,1])#peak at 13.2
x_182_peak_eV,y_182_peak=np.loadtxt('data/Y24906L.txt', unpack=True, usecols= [0,1])#peak at 18.2

x_ex0 = x#np.arange(0, 10, 0.4)
y_ex0 = y#,y_042_peak,y_062_peak,y_082_peak,y_102_peak,y_132_peak,y_182_peak]
z_ex0 = 0#[0.0, 4.2, 6.2, 8.2, 10.2, 13.2, 18.2]
vert_ex0 = []
#for z in zs:
#ys[0], ys[-1] = 0, 0
vert_ex0.append(list(zip(x_ex0, y_ex0)))
poly_ex0 = PolyCollection(vert_ex0, facecolors = cc('#FFFFFF'))#, cc('c'), cc('r'), cc('g'), cc('y'), cc('m'), cc('k')])
poly_ex0.set_alpha(0.9)
ax.add_collection3d(poly_ex0, zs=z_ex0, zdir='y')

x_0 = x_nopeak_eV#np.arange(0, 10, 0.4)
y_0 = y_nopeak#,y_042_peak,y_062_peak,y_082_peak,y_102_peak,y_132_peak,y_182_peak]
z_0 = 10.0#[0.0, 4.2, 6.2, 8.2, 10.2, 13.2, 18.2]
vert_0 = []
#for z in zs:
#ys[0], ys[-1] = 0, 0
vert_0.append(list(zip(x_0, y_0)))
poly_0 = PolyCollection(vert_0,  facecolors = cc('#FFDD00'))#, cc('c'), cc('r'), cc('g'), cc('y'), cc('m'), cc('k')])
poly_0.set_alpha(0.9)
ax.add_collection3d(poly_0, zs=z_0, zdir='y')

x_1 = x_nopeak_eV#np.arange(0, 10, 0.4)
y_1 = y_042_peak#,y_042_peak,y_062_peak,y_082_peak,y_102_peak,y_132_peak,y_182_peak]
z_1 = 20.0#[0.0, 4.2, 6.2, 8.2, 10.2, 13.2, 18.2]
vert_1 = []
#for z in zs:
#ys[0], ys[-1] = 0, 0
vert_1.append(list(zip(x_1, y_1)))
poly_1 = PolyCollection(vert_1,  facecolors = cc('#FF9900'))#[0.78,0.28,0.8]))#, cc('c'), cc('r'), cc('g'), cc('y'), cc('m'), cc('k')])
poly_1.set_alpha(0.9)
ax.add_collection3d(poly_1, zs=z_1, zdir='y')

x_2 = x_nopeak_eV#np.arange(0, 10, 0.4)
y_2 = y_062_peak#,y_042_peak,y_062_peak,y_082_peak,y_102_peak,y_132_peak,y_182_peak]
z_2 = 30.0#[0.0, 4.2, 6.2, 8.2, 10.2, 13.2, 18.2]
vert_2 = []
#for z in zs:
#ys[0], ys[-1] = 0, 0
vert_2.append(list(zip(x_2, y_2)))
poly_2 = PolyCollection(vert_2, facecolors = cc('#FF5500'))#[0.5,0,0]))#, cc('c'), cc('r'), cc('g'), cc('y'), cc('m'), cc('k')])
poly_2.set_alpha(0.9)
ax.add_collection3d(poly_2, zs=z_2, zdir='y')

x_3 = x_nopeak_eV#np.arange(0, 10, 0.4)
y_3 = y_082_peak#,y_042_peak,y_062_peak,y_082_peak,y_102_peak,y_132_peak,y_182_peak]
z_3 = 40.0#[0.0, 4.2, 6.2, 8.2, 10.2, 13.2, 18.2]
vert_3 = []
#for z in zs:
#ys[0], ys[-1] = 0, 0
vert_3.append(list(zip(x_3, y_3)))
poly_3 = PolyCollection(vert_3, facecolors = cc('#FF0000'))#[0.53,0.94,0]))#, cc('c'), cc('r'), cc('g'), cc('y'), cc('m'), cc('k')])
poly_3.set_alpha(0.9)
ax.add_collection3d(poly_3, zs=z_3, zdir='y')

x_4 = x_nopeak_eV#np.arange(0, 10, 0.4)
y_4 = y_102_peak#,y_042_peak,y_062_peak,y_082_peak,y_102_peak,y_132_peak,y_182_peak]
z_4 = 50.0#[0.0, 4.2, 6.2, 8.2, 10.2, 13.2, 18.2]
vert_4 = []
#for z in zs:
#ys[0], ys[-1] = 0, 0
vert_4.append(list(zip(x_4, y_4)))
poly_4 = PolyCollection(vert_4,facecolors = cc('#663300'))#[0.255,0.149,0.014]))#, cc('c'), cc('r'), cc('g'), cc('y'), cc('m'), cc('k')])
poly_4.set_alpha(0.9)
ax.add_collection3d(poly_4, zs=z_4, zdir='y')

#x_6 = x_nopeak_eV#np.arange(0, 10, 0.4)
#y_6 = y_182_peak#,y_042_peak,y_062_peak,y_082_peak,y_102_peak,y_132_peak,y_182_peak]
#z_6 = 18.2#[0.0, 4.2, 6.2, 8.2, 10.2, 13.2, 18.2]
#vert_6 = []
##for z in zs:
##ys[0], ys[-1] = 0, 0
#vert_6.append(list(zip(x_6, y_6)))
#poly_6 = PolyCollection(vert_6, facecolors = cc('#990000'))#'k'))#, cc('c'), cc('r'), cc('g'), cc('y'), cc('m'), cc('k')])
#poly_6.set_alpha(0.9)
#ax.add_collection3d(poly_6, zs=z_6, zdir='y')

x_5 = x_nopeak_eV#np.arange(0, 10, 0.4)
y_5 = y_132_peak#,y_042_peak,y_062_peak,y_082_peak,y_102_peak,y_132_peak,y_182_peak]
z_5 = 60.0#[0.0, 4.2, 6.2, 8.2, 10.2, 13.2, 18.2]
vert_5 = []
#for z in zs:
#ys[0], ys[-1] = 0, 0
vert_5.append(list(zip(x_5, y_5)))
poly_5 = PolyCollection(vert_5, facecolors = cc('#660000'))#[0,0.4,0.4]))#, cc('c'), cc('r'), cc('g'), cc('y'), cc('m'), cc('k')])
poly_5.set_alpha(0.9)
ax.add_collection3d(poly_5, zs=z_5, zdir='y')

#ax.set_xlabel(r'$\hbar\omega\,\,[eV]$')
ax.set_xlim3d(0, 50)
#ax.set_ylabel('Method')
ax.set_ylim3d(60, 0)
#ax.set_zlabel(r'$\cal{Im}$(Dielectric Response)')
ax.set_zlim3d(0, 4)
plt.show()
