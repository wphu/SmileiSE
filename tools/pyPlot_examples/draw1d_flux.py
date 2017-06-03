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


	
ylabel_x = -0.058
	
	
##inite the fig of matplotlib
fig=plt.figure(figsize=(10,8))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.6,hspace=0.55)


t = 14
dt = 1.0e-12 * 1000000 / 1.0e-6   # unit is us

##read data from file
f=h5.File("data_global.h5")
print f.keys()

group = f['/Fields']


val1 = f["/Diagnostic/particleNumber"]
val1 = val1[...]

nx = (val1.shape)[0]
x = np.linspace(0,nx*dt,nx)
print nx


#============Total particle number=======================================

val1 = f["/Diagnostic/particleNumber"]
val1 = val1[...]
val1_1d = np.transpose(val1[:, 0, 0, 0])
#val1_1d.astype(float)

val2 = f["/Diagnostic/particleNumber"]
val2 = val2[...]
val2_1d = np.transpose(val2[:, 0, 0, 1])
#val2_1d.astype(float)

sp_temp1=fig.add_subplot(3,1,1)


sp_temp1.yaxis.set_major_formatter(yformatter)

cf_temp1=sp_temp1.plot(x, val1_1d, label='Electron')
cf_temp1=sp_temp1.plot(x, val2_1d, label=r'$D^+ ion$')

xmin = 0
xmax = x.max()
ymin = 0.0
ymax = 1.2 * max( val1_1d.max(), val2_1d.max() )

sp_temp1.set_xlim((xmin, xmax))
sp_temp1.set_ylim((ymin, ymax))
sp_temp1.set_ylabel(r'Macro particle number')
sp_temp1.yaxis.set_label_coords(ylabel_x, 0.5)

#sp_temp1.ticklabel_format(style='sci')
sp_temp1.grid(True)
sp_temp1.legend(loc = 1)

sp_temp1.annotate('(a)', xy=get_axis_limits(sp_temp1), annotation_clip=False)

##============Particle flux======================================================
val1 = f["/Diagnostic/particleFlux"]
val1 = val1[...]
val1_1d = np.transpose(val1[:, 0, 0, 0])

val2 = f["/Diagnostic/particleFlux"]
val2 = val2[...]
val2_1d = np.transpose(val2[:, 0, 1, 0])


sp_temp1=fig.add_subplot(3,1,2)

sp_temp1.yaxis.set_major_formatter(yformatter)

cf_temp1=sp_temp1.plot(x, val1_1d, label='Electron')
cf_temp1=sp_temp1.plot(x, val2_1d, label=r'$D^+ ion$')


ymin = 0.0
ymax = 1.2 * max( val1_1d.max(), val2_1d.max() )

sp_temp1.set_xlim((xmin, xmax))
sp_temp1.set_ylim((ymin, ymax))
sp_temp1.set_ylabel(r'Particle flux $\mathit{(m^{-2}s^{-1})}$')
sp_temp1.yaxis.set_label_coords(ylabel_x, 0.5)

sp_temp1.grid(True)
sp_temp1.legend(loc = 1)

sp_temp1.annotate('(b)', xy=get_axis_limits(sp_temp1), annotation_clip=False)

##============Eenergy flux======================================================
val1 = f["/Diagnostic/heatFlux"]
val1 = val1[...]
val1_1d = np.transpose(val1[:, 0, 0, 0])

val2 = f["/Diagnostic/heatFlux"]
val2 = val2[...]
val2_1d = np.transpose(val2[:, 0, 1, 0])


sp_temp1=fig.add_subplot(3,1,3)

sp_temp1.yaxis.set_major_formatter(yformatter)

cf_temp1=sp_temp1.plot(x, val1_1d, label='Electron')
cf_temp1=sp_temp1.plot(x, val2_1d, label=r'$D^+ ion$')


ymin = 0.0
ymax = 1.2 * max( val1_1d.max(), val2_1d.max() )

sp_temp1.set_xlim((xmin, xmax))
sp_temp1.set_ylim((ymin, ymax))
sp_temp1.set_xlabel(r'time $(\mu s)$')
sp_temp1.set_ylabel(r'Heat flux $(Wm^{-2})$')
sp_temp1.yaxis.set_label_coords(ylabel_x, 0.5)
sp_temp1.grid(True)
sp_temp1.legend(loc = 1)

sp_temp1.annotate('(c)', xy=get_axis_limits(sp_temp1), annotation_clip=False)

fig.savefig("flux.pdf")
#fig.show()       #when the program finishes,the figure disappears
#plt.axis('equal')
#plt.show()         #The command is OK


##1d plot=============================
##ax1=fig.add_subplot(2,1,1)
#flux=f["/1d/pflux"]
#flux=flux[...]

#line1,=ax1.plot(flux[0,:],flux[1,:])
#line2,=ax1.plot(flux[0,:],flux[2,:])
