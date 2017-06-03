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


mpl.rcParams['font.size'] = 16
mpl.rcParams['axes.linewidth'] = 2.0
#mpl.rcParams['font.weight'] = 'bold'

mpl.rcParams['xtick.major.size'] = 2
mpl.rcParams['ytick.major.size'] = 2

mpl.rcParams['lines.linewidth'] = 2.0

#mpl.rcParams['grid.linestyle'] = ":"
mpl.rcParams['grid.linestyle'] = ":"
mpl.rcParams['grid.color'] = "black"

def get_axis_limits(ax, x_scale=0, y_scale=1.18):
    return ax.get_xlim()[1]*x_scale, ax.get_ylim()[1]*y_scale



##inite the fig of matplotlib
fig=plt.figure(figsize=(10,8))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.55)

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
xmax = x.max() * 0.01/4
x_sheath = 0.08



##============potential======================================================
sp_temp1=fig.add_subplot(3,1,1)

val = f["/Fields/Phi_global_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])

cf_temp1=sp_temp1.plot(x, val_1d, label = "Electric potential $(V)$", color='#1f77b4')

sp_temp1.set_ylabel("Electric potential $(V)$", color='#1f77b4')
sp_temp1.tick_params('y', colors='#1f77b4')



#double y axis
sp_temp2 = sp_temp1.twinx()
val = f["/Fields/Ex_global_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])

sp_temp2.yaxis.set_major_formatter(yformatter)
cf_temp1=sp_temp2.plot(x, val_1d, label = "Electric field (x)", color='#ff7f0e')

sp_temp2.set_ylabel("$E_x \ (V/m)$", color='#ff7f0e')
sp_temp2.tick_params('y', colors='#ff7f0e')

#lines1, labels1 = sp_temp1.get_legend_handles_labels()
#lines2, labels2 = sp_temp2.get_legend_handles_labels()
#sp_temp2.legend(lines1 + lines2, labels1 + labels2, loc = 1, fancybox = False)

sp_temp1.grid(True)


#sp_temp1.axis([x.min(),x.max(),val_1d.min(),val_1d.max()])
#sp_temp1.set_yticks(np.arange(0,y.max(),100))
sp_temp1.set_xlim((xmin, xmax))
#sp_temp1.set_xlabel('x(mm)')
#sp_temp1.set_ylabel('Electric potential (V)')
#sp_temp2.set_ylabel('Electric field (V/m)')

sp_temp1.annotate('(a)', xy=get_axis_limits(sp_temp1), annotation_clip=False)

##============rho======================================================
sp_temp1=fig.add_subplot(3,1,2)

val = f["/Fields/Rho_global_e_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])

sp_temp1.yaxis.set_major_formatter(yformatter)
cf_temp1=sp_temp1.plot(x, val_1d, label = "Electron")

ymin = 0.0
ymax = val_1d.max()

val = f["/Fields/Rho_global_D1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = r'$D^+ ion$')

sp_temp1.grid(True)
sp_temp1.legend(loc = 1)
sp_temp1.set_xlim((xmin, xmax))
sp_temp1.set_ylim((ymin, ymax))
#sp_temp1.set_yticks(np.arange(0,y.max(),100))
#sp_temp1.set_xlabel('x(mm)')
sp_temp1.set_ylabel('Number density $(m^{-3})$')

sp_temp1.annotate('(b)', xy=get_axis_limits(sp_temp1), annotation_clip=False)

##============ Temperature ======================================================
sp_temp1=fig.add_subplot(3,1,3)

val = f["/Fields/T_global_e_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = "Electron")

val = f["/Fields/T_global_D1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = r'$D^+ ion$')

ymin = 0.0
ymax = val_1d.max()
# Plot a line
sp_temp1.axvline(x = x_sheath, ymin = 0.0, ymax = 3.8, c="red",zorder=0, clip_on=False)


sp_temp1.grid(True)
sp_temp1.legend(loc = 1)
sp_temp1.set_xlim((xmin, xmax))
sp_temp1.set_ylim((ymin, ymax))
#sp_temp1.set_yticks(np.arange(0,y.max(),100))
sp_temp1.set_xlabel('x $(mm)$')
sp_temp1.set_ylabel('Temperature $(eV)$')

sp_temp1.annotate('(c)', xy=get_axis_limits(sp_temp1), annotation_clip=False)

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



# ion sound speed
Va = math.sqrt( 60.0 * 1.602e-19 / (2.0 * 1.67262158e-27) ) / 1.0e6
# Plot ion sound speed
#sp_temp1.axhline(y = Va, xmin = 0.0, xmax = 1.0, c="red",zorder=0, clip_on=False)
sp_temp1.axhline(y = -Va, xmin = 0.0, xmax = 1.0, c="red",zorder=0, clip_on=False)

sp_temp1.grid(True)
sp_temp1.legend()
sp_temp1.set_xlim((xmin, xmax))
#sp_temp1.set_yticks(np.arange(0,y.max(),100))
sp_temp1.set_xlabel('x(mm)')
sp_temp1.set_ylabel('Vp (1.0e6 m/s)')
'''


fig.savefig("Sheath.pdf", dpi = 300)
##fig.show()       #when the program finishes,the figure disappears
#plt.axis('equal')
plt.show()         #The command is OK


##1d plot=============================
##ax1=fig.add_subplot(2,1,1)
#flux=f["/1d/pflux"]
#flux=flux[...]

#line1,=ax1.plot(flux[0,:],flux[1,:])
#line2,=ax1.plot(flux[0,:],flux[2,:])
