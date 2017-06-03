##>>>The code is used to read data from hdf5 file
##>>>and plot on the screen and output figure file using matplotlib-python

import Tkinter as tk
from Tkinter import *

import matplotlib
matplotlib.use('Agg')


from matplotlib.ticker import MaxNLocator
from matplotlib import cm

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure

from matplotlib.ticker import ScalarFormatter
yformatter = ScalarFormatter()
yformatter.set_powerlimits((-3,3))

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker


from numpy import arange, sin, pi
import numpy as np

import h5py as h5
import math



font={	'family' : 'sans-serif',
	'weight' : 'bold',
	'size' : 8,
	}

#mpl.rcParams['text.usetex'] = True
#mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['font.family'] = 'sans-serif'
#mpl.rcParams['mathtext.default'] = 'regular'
#mpl.rcParams['mathtext.default'] = 'it'
mpl.rcParams['mathtext.fontset'] = 'stix'
#mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.it'] = 'serif'
#mpl.rcParams['pdf.fonttype'] = 3


mpl.rcParams['font.size'] = 18
mpl.rcParams['axes.linewidth'] = 2.0
#mpl.rcParams['font.weight'] = 'bold'

mpl.rcParams['xtick.major.size'] = 2
mpl.rcParams['ytick.major.size'] = 2

mpl.rcParams['lines.linewidth'] = 2.0

#mpl.rcParams['grid.linestyle'] = ":"
mpl.rcParams['grid.linestyle'] = ":"
mpl.rcParams['grid.color'] = "black"

def get_axis_limits(ax, x_scale=0, y_scale=1.16):
    return ax.get_xlim()[1]*x_scale, ax.get_ylim()[1]*y_scale



##inite the fig of matplotlib
fig=plt.figure(figsize=(10,8))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.5)

t = 19


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
xmax = x.max() * 0.1
ymin = 0.0


##============rho======================================================
sp_temp1=fig.add_subplot(2,1,1)
sp_temp1.yaxis.set_major_formatter(yformatter)

val = f["/Fields/Rho_global_e_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = r'$Electron$')

val = f["/Fields/Rho_global_D1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = r'$D^+ \ ion$')


val = f["/Fields/Rho_global_C_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = r'$C \ atom$')



val = f["/Fields/Rho_global_C1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = r'$C^+ \ ion$')

val = f["/Fields/Rho_global_C2_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = r'$C^{+2} \ ion$')

val = f["/Fields/Rho_global_C3_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = r'$C^{+3} \ ion$')


sp_temp1.grid(True)
sp_temp1.legend(loc = 1, framealpha = 1.0)
sp_temp1.set_xlim((xmin, xmax))
sp_temp1.set_ylim((ymin))
#sp_temp1.set_yticks(np.arange(0,y.max(),100))
#sp_temp1.set_xlabel('x(mm)')
sp_temp1.set_ylabel(r'Number density $(m^{-3})$')

sp_temp1.annotate('(a)', xy=get_axis_limits(sp_temp1), annotation_clip=False)

##============ Temperature ======================================================
sp_temp1=fig.add_subplot(2,1,2)

val = f["/Fields/T_global_e_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = r'Electron')

val = f["/Fields/T_global_D1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = r'$D^+ \ ion$')
ymax = val_1d.max() * 1.1


val = f["/Fields/T_global_C_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = r'C atom')


val = f["/Fields/T_global_C1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = r'$C^+ \ ion$')


val = f["/Fields/T_global_C2_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = r'$C^{+2} \ ion$')


val = f["/Fields/T_global_C3_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = r'$C^{+3} \ ion$')

ymin = 0

sp_temp1.grid(True)
#sp_temp1.legend(loc = 1)
sp_temp1.set_xlim((xmin, xmax))
sp_temp1.set_ylim((ymin, ymax))
#sp_temp1.set_yticks(np.arange(0,y.max(),100))
sp_temp1.set_xlabel(r'x $(mm)$')
sp_temp1.set_ylabel(r'Temperature $(eV)$')

sp_temp1.annotate('(b)', xy=get_axis_limits(sp_temp1), annotation_clip=False)

fig.savefig("nT_C.pdf", dpi = 300)
##fig.show()       #when the program finishes,the figure disappears
#plt.axis('equal')
#plt.show()         #The command is OK


##1d plot=============================
##ax1=fig.add_subplot(2,1,1)
#flux=f["/1d/pflux"]
#flux=flux[...]

#line1,=ax1.plot(flux[0,:],flux[1,:])
#line2,=ax1.plot(flux[0,:],flux[2,:])
