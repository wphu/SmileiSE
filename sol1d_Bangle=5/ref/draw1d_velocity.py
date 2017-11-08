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

def get_axis_limits(ax, x_scale=-0.095, y_scale=1.03):
    return ax.get_xlim()[1]*x_scale, ax.get_ylim()[1] + ( ax.get_ylim()[1] - ax.get_ylim()[0] ) * (y_scale - 1.0)



##inite the fig of matplotlib
fig=plt.figure(figsize=(10,8))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.4)

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
xmax = x.max()

# ion sound speed
Va1 = math.sqrt( 90.3 * 1.602e-19 / (2.0 * 1.67262158e-27) ) 
Va2 = math.sqrt( 120.0 * 1.602e-19 / (2.0 * 1.67262158e-27) )

##============ Parallel Velocity ======================================================
sp_temp1=fig.add_subplot(2,1,1)
sp_temp1.yaxis.set_major_formatter(yformatter)

val = f["/Fields/Vparallel_global_D1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
val1_1d = val_1d
cf_temp1=sp_temp1.plot(x, val1_1d, label = r"$\mathrm{D^+}$ ion")

val = f["/Fields/Vparallel_global_e_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
val0_1d = val_1d
cf_temp1=sp_temp1.plot(x, val0_1d, label = r"Electron")


'''
x0_ax_inset0 = 0.0
y0_ax_inset0 = val0_1d.max() * 0.1
width0_ax_inset0 = 8.0
height0_ax_inset0 = val0_1d.max() * 0.8
ax_inset0 = inset_axes(sp_temp1, width = "60%", height = "80%" , bbox_to_anchor = (x0_ax_inset0, y0_ax_inset0, width0_ax_inset0, height0_ax_inset0), bbox_transform = sp_temp1.transData)

ax_inset0.yaxis.set_major_formatter(yformatter)
ax_inset0.plot(x, val1_1d)
ax_inset0.plot(x, val0_1d)
ax_inset0.set_xlim(0.0, 0.2)
#ax_inset0.set_ylim(0.0, 0.3)
'''

# Plot ion sound speed
#sp_temp1.axhline(y = Va1, xmin = 0.0, xmax = 1.0, c="red",zorder=0, clip_on=False)
#sp_temp1.axhline(y = -Va1, xmin = 0.0, xmax = 1.0, c="red",zorder=0, clip_on=False)
#sp_temp1.axhline(y = Va2, xmin = 0.0, xmax = 1.0, c="blue",zorder=0, clip_on=False)
#sp_temp1.axhline(y = -Va2, xmin = 0.0, xmax = 1.0, c="blue",zorder=0, clip_on=False)

sp_temp1.legend()
xmin = 0.0
xmax = 0.25
sp_temp1.set_xlim((xmin, xmax))
sp_temp1.set_xlim((xmin, xmax))
#sp_temp1.set_yticks(np.arange(0,y.max(),100))
sp_temp1.set_xlabel(r"$x\ \mathrm{(mm)}$", fontsize = label_fontsize)
sp_temp1.set_ylabel(r"$V\ \mathrm{(m/s)}$", fontsize = label_fontsize)

sp_temp1.annotate(r"$\mathbf{(a)}$", xy=get_axis_limits(sp_temp1), annotation_clip=False)

##============ Ion Parallel Velocity ======================================================
sp_temp1=fig.add_subplot(2,1,2)
sp_temp1.yaxis.set_major_formatter(yformatter)


val = f["/Fields/Vparallel_global_D1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
val1_1d = val_1d
cf_temp1=sp_temp1.plot(x, val1_1d, label = "D+1")



# Plot ion sound speed
#sp_temp1.axhline(y = Va1, xmin = 0.0, xmax = 1.0, c="red",zorder=0, clip_on=False)
sp_temp1.axhline(y = -Va1, xmin = 0.0, xmax = 1.0, c="red",zorder=0, clip_on=False)
#sp_temp1.axhline(y = Va2, xmin = 0.0, xmax = 1.0, c="blue",zorder=0, clip_on=False)
#sp_temp1.axhline(y = -Va2, xmin = 0.0, xmax = 1.0, c="blue",zorder=0, clip_on=False)

#sp_temp1.legend()
xmin = 0.0
xmax = 0.25
ymin = val1_1d.min() * 1.1
ymax = 0.0
sp_temp1.set_xlim((xmin, xmax))
sp_temp1.set_ylim((ymin, ymax))
#sp_temp1.set_yticks(np.arange(0,y.max(),100))
sp_temp1.set_xlabel(r"$x\ \mathrm{(mm)}$", fontsize = label_fontsize)
sp_temp1.set_ylabel(r"$V_{\mathrm{D^+}}\ \mathrm{(m/s)}$", fontsize = label_fontsize)

major_ticks = np.arange(0, -1.2e5, -0.5e5)                                                                                      
sp_temp1.set_yticks(major_ticks)   

sp_temp1.annotate(r"$\mathbf{(b)}$", xy=get_axis_limits(sp_temp1), annotation_clip=False)


'''
##============ x direction Velocity ======================================================
sp_temp1=fig.add_subplot(2,1,2)

val = f["/Fields/Vx_global_e_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d / 1.0e6
cf_temp1=sp_temp1.plot(x, val_1d, label = "Electron")

val = f["/Fields/Vy_global_D1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d / 1.0e6
cf_temp1=sp_temp1.plot(x, val_1d, label = "D+1")

sp_temp1.legend()
sp_temp1.set_xlim((xmin, xmax))
#sp_temp1.set_yticks(np.arange(0,y.max(),100))
sp_temp1.set_xlabel('x(mm)')
sp_temp1.set_ylabel('Vx (1.0e6 m/s)')
'''


fig.savefig("velocity.pdf", dpi = 300)
##fig.show()       #when the program finishes,the figure disappears
#plt.axis('equal')
plt.show()         #The command is OK


##1d plot=============================
##ax1=fig.add_subplot(2,1,1)
#flux=f["/1d/pflux"]
#flux=flux[...]

#line1,=ax1.plot(flux[0,:],flux[1,:])
#line2,=ax1.plot(flux[0,:],flux[2,:])
