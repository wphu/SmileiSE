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



font={	'family' : 'serif',
	'weight' : 'bold',
	'size' : 8,
	}

#mpl.rcParams['text.usetex'] = True
#mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['font.family'] = 'serif'
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

mpl.rcParams['lines.linewidth'] = 3.0

#mpl.rcParams['grid.linestyle'] = ":"
mpl.rcParams['grid.linestyle'] = ":"
mpl.rcParams['grid.color'] = "black"

label_fontsize = 21
legend_fontsize = 15


def get_axis_limits(ax, x_scale=-0.1, y_scale=1.03):
    return ax.get_xlim()[0] + (ax.get_xlim()[1] - ax.get_xlim()[0]) * x_scale, ax.get_ylim()[1] + ( ax.get_ylim()[1] - ax.get_ylim()[0] ) * (y_scale - 1.0)



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


dx = 0.5e-5  # unit (m)
x = np.linspace(0, nx * dx, nx)

amplification_factor = 80.0
x = x * amplification_factor

xmin = x.min()
xmax = x.max() * 0.2
ymin = 0.0


##============rho======================================================
sp_temp1=fig.add_subplot(2,1,1)
sp_temp1.yaxis.set_major_formatter(yformatter)

val = f["/Fields/Rho_global_e_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = r'Electron')

val = f["/Fields/Rho_global_D1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = r'$\mathrm{D^+ \ ion}$')


val = f["/Fields/Rho_global_C_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = r'$\mathrm{C \ atom}$')



val = f["/Fields/Rho_global_C1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = r'$\mathrm{C^+ \ ion}$')

val = f["/Fields/Rho_global_C2_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = r'$\mathrm{C^{+2} \ ion}$')

val = f["/Fields/Rho_global_C3_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = r'$\mathrm{C^{+3} \ ion}$')


sp_temp1.grid(True)
sp_temp1.legend(loc = 1, framealpha = 1.0, fontsize = legend_fontsize)
sp_temp1.set_xlim((xmin, xmax))
sp_temp1.set_ylim((ymin))

major_ticks = np.arange(0, 2.1e19, 0.5e19)                                              
#minor_ticks = np.arange(0, 31, 5)                                          
sp_temp1.set_yticks(major_ticks)                                                       
#sp_temp1.set_yticks(minor_ticks, minor=True)  

#sp_temp1.set_yticks(np.arange(0,y.max(),100))
#sp_temp1.set_xlabel('x(mm)')
sp_temp1.set_ylabel(r"$n\ \mathrm{(m^{-3})}$", fontsize = label_fontsize)

sp_temp1.annotate(r"$\mathbf{(a)}$", xy=get_axis_limits(sp_temp1), annotation_clip=False)

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
sp_temp1.set_ylim((ymin, 66))

major_ticks = np.arange(0, 71.0, 20.0)                                              
#minor_ticks = np.arange(0, 31, 5)                                          
sp_temp1.set_yticks(major_ticks)                                                       
#sp_temp1.set_yticks(minor_ticks, minor=True)  


#sp_temp1.set_yticks(np.arange(0,y.max(),100))
sp_temp1.set_xlabel(r"$x\ \mathrm{(m)}$", fontsize = label_fontsize)
sp_temp1.set_ylabel(r"$T\ \mathrm{(eV)}$", fontsize = label_fontsize)

sp_temp1.annotate(r"$\mathbf{(b)}$", xy=get_axis_limits(sp_temp1), annotation_clip=False)

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
