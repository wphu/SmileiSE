##>>>The code is used to read data from hdf5 file
##>>>and plot on the screen and output figure file using matplotlib-python

import Tkinter as tk
from Tkinter import *

import matplotlib
matplotlib.use('Agg')


import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
from matplotlib import cm


from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure

import numpy as np
from numpy import arange, sin, pi

import ConfigParser

import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.ticker import ScalarFormatter

yformatter = ScalarFormatter()
yformatter.set_powerlimits((-3,3))




font={	'family' : 'sans-serif',
	'weight' : 'bold',
	'size' : 8,
	}

mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['mathtext.default'] = 'regular'

mpl.rcParams['font.size'] = 16
mpl.rcParams['axes.linewidth'] = 2.0
#mpl.rcParams['font.weight'] = 'bold'

mpl.rcParams['xtick.major.size'] = 2
mpl.rcParams['ytick.major.size'] = 2

mpl.rcParams['lines.linewidth'] = 2.0

mpl.rcParams['grid.linestyle'] = "--"
mpl.rcParams['grid.color'] = "black"




##inite the fig of matplotlib
fig=plt.figure(figsize=(10,8))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.4)
t = 15


##read data from file
f=h5.File("data_global.h5")
print f.keys()

group = f['/Fields']
dims = group.attrs['dims_global']
dims = dims[...]

print dims

nx = dims[3]


dx = 0.5e-2  # unit (mm)
x = np.linspace(0, nx * dx, nx)

xmin = x.min()
xmax = x.max()

x1=5000
x2=4500
x3=12

# ion sound speed
Va = - math.sqrt( 25.0 * 1.602e-19 / (2.0 * 1.67262158e-27) ) / 1.0e6


##============potential======================================================
sp_temp1=fig.add_subplot(3,1,1)

val = f["/Fields/Phi_global_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = "Electric potential")
print "potential max: ", val_1d.max()
print "potential: ",val_1d[x1], val_1d[x2], val_1d[x3]

#double y axis
sp_temp2 = sp_temp1.twinx()
val = f["/Fields/Ex_global_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])

sp_temp2.yaxis.set_major_formatter(yformatter)
cf_temp1=sp_temp2.plot(x, val_1d, label = "Electric field (x)", color = "green")

sp_temp1.grid(True)

sp_temp1.legend()
sp_temp2.legend()

legend1=sp_temp1.legend(loc=(.6,.76),fontsize=16)
legend2=sp_temp2.legend(loc=(.6,.56),fontsize=16)
#sp_temp1.axis([x.min(),x.max(),val_1d.min(),val_1d.max()])
#sp_temp1.set_yticks(np.arange(0,y.max(),100))
sp_temp1.set_xlim((xmin, xmax))
sp_temp1.set_xlabel('x(mm)')
sp_temp1.set_ylabel('Electric potential (V)')
sp_temp2.set_ylabel('Electric field (V/m)')


##============rho======================================================
sp_temp1=fig.add_subplot(3,1,2)

val = f["/Fields/Rho_global_e_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = "Electron")
print "Electron density: ",val_1d[x1], val_1d[x2], val_1d[x3]



val = f["/Fields/Rho_global_D1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = "D+1")

sp_temp1.grid(True)

sp_temp1.legend()
sp_temp1.set_xlim((xmin, xmax))
#sp_temp1.set_yticks(np.arange(0,y.max(),100))
sp_temp1.set_xlabel('x(mm)')
sp_temp1.set_ylabel('Number density')



##============ Temperature ======================================================
sp_temp1=fig.add_subplot(3,1,3)

val = f["/Fields/T_global_e_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = "Electron")
print "Electron Temperature: ",val_1d[x1], val_1d[x2], val_1d[x3]

sp_temp1.grid(True)

val = f["/Fields/T_global_D1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = "D+1")
print "D+1 ion temperature: ",val_1d[x1], val_1d[x2], val_1d[x3]

sp_temp1.legend()


ymin = 0
ymax = val_1d.max() * 1.2
sp_temp1.set_xlim((xmin, xmax))
sp_temp1.set_ylim((ymin, ymax))
#sp_temp1.set_yticks(np.arange(0,y.max(),100))
sp_temp1.set_xlabel('x(mm)')
sp_temp1.set_ylabel('Temperature (eV)')


##============ Velocity ======================================================
'''
sp_temp1=fig.add_subplot(4,1,4)

val = f["/Fields/Vparallel_global_e_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d / 1.0e6
cf_temp1=sp_temp1.plot(x, val_1d, label = "Electron")


val = f["/Fields/Vparallel_global_D1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d / 1.0e6
cf_temp1=sp_temp1.plot(x, val_1d, label = "D+1")

# Plot ion sound speed
sp_temp1.axhline(y = Va, xmin = 0.0, xmax = 1.0, c="red",zorder=0, clip_on=False)
sp_temp1.axhline(y = -Va, xmin = 0.0, xmax = 1.0, c="red",zorder=0, clip_on=False)

sp_temp1.legend()
sp_temp1.set_xlim((xmin, xmax))
#sp_temp1.set_yticks(np.arange(0,y.max(),100))
sp_temp1.set_xlabel('x(mm)')
sp_temp1.set_ylabel('Vp (1.0e6 m/s)')
'''



fig.savefig("Profiles.png", dpi = 300)
##fig.show()       #when the program finishes,the figure disappears
#plt.axis('equal')
#plt.show()         #The command is OK


##1d plot=============================
##ax1=fig.add_subplot(2,1,1)
#flux=f["/1d/pflux"]
#flux=flux[...]

#line1,=ax1.plot(flux[0,:],flux[1,:])
#line2,=ax1.plot(flux[0,:],flux[2,:])
