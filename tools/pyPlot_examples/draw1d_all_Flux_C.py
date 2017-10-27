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

label_fontsize = 23
legend_fontsize = 23


def get_axis_limits(ax, x_scale=-0.1, y_scale=1.03):
    return ax.get_xlim()[0] + (ax.get_xlim()[1] - ax.get_xlim()[0]) * x_scale, ax.get_ylim()[1] + ( ax.get_ylim()[1] - ax.get_ylim()[0] ) * (y_scale - 1.0)



##inite the fig of matplotlib
fig=plt.figure(figsize=(10,8))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.55)

t = 19




##read data from file
f=h5.File("ref/data_global.h5")
print f.keys()

group = f['/Fields']
dims = group.attrs['dims_global']
dims = dims[...]

print dims

nx = dims[3]


dx = 0.5e-2  # unit (mm)
x = np.zeros(4)
y1 = np.zeros(4)
y2 = np.zeros(4)
y3 = np.zeros(4)
pFlux = np.zeros(4)		# particle flux
hFlux = np.zeros(4)		# heat flux

x[0] = 0.0
x[1] = 1.0
x[2] = 2.0
x[3] = 3.0

x = x * 1.0e22

#xmin = x.min()
#xmax = x.max()
ymin = 0.0
ymax = 6.4e6



sp_temp1=fig.add_subplot(1,1,1)
sp_temp1.yaxis.set_major_formatter(yformatter)

# ================================= heat flux =========================================
##============Ref =======
f=h5.File("ref/data_global.h5")

val1 = f["/Diagnostic/heatFlux"]
val1 = val1[...]
val1_1d = np.transpose(val1[:, 0, 0, 0])
y1[0] = val1_1d[t]

val2 = f["/Diagnostic/heatFlux"]
val2 = val2[...]
val2_1d = np.transpose(val2[:, 0, 1, 0])
y2[0] = val2_1d[t]

#C1 ion
y3[0] = 0.0

##============ Case1 =========

f=h5.File("IC1.0/data_global.h5")

val1 = f["/Diagnostic/heatFlux"]
val1 = val1[...]
val1_1d = np.transpose(val1[:, 0, 0, 0])
y1[1] = val1_1d[t]

val2 = f["/Diagnostic/heatFlux"]
val2 = val2[...]
val2_1d = np.transpose(val2[:, 0, 1, 0])
y2[1] = val2_1d[t]

val2 = f["/Diagnostic/heatFlux"]
val2 = val2[...]
val2_1d = np.transpose(val2[:, 0, 4, 0])
y3[1] = val2_1d[t]

##============ Case2 ============
f=h5.File("IC2.0/data_global.h5")

val1 = f["/Diagnostic/heatFlux"]
val1 = val1[...]
val1_1d = np.transpose(val1[:, 0, 0, 0])
y1[2] = val1_1d[t]

val2 = f["/Diagnostic/heatFlux"]
val2 = val2[...]
val2_1d = np.transpose(val2[:, 0, 1, 0])
y2[2] = val2_1d[t]

val2 = f["/Diagnostic/heatFlux"]
val2 = val2[...]
val2_1d = np.transpose(val2[:, 0, 4, 0])
y3[2] = val2_1d[t]

##============ Case3 ===========
f=h5.File("IC3.0/data_global.h5")

val1 = f["/Diagnostic/heatFlux"]
val1 = val1[...]
val1_1d = np.transpose(val1[:, 0, 0, 0])
y1[3] = val1_1d[t]

val2 = f["/Diagnostic/heatFlux"]
val2 = val2[...]
val2_1d = np.transpose(val2[:, 0, 1, 0])
y2[3] = val2_1d[t]

val2 = f["/Diagnostic/heatFlux"]
val2 = val2[...]
val2_1d = np.transpose(val2[:, 0, 4, 0])
y3[3] = val2_1d[t]


hFlux = y1 + y2 + y3
hFlux_norm  = hFlux[0]	# used to normalize hFlux
print "hFlux: ", hFlux_norm
hFlux = hFlux / hFlux_norm




# ========================================== particle flux =======================
##============Ref =============
f=h5.File("ref/data_global.h5")

val1 = f["/Diagnostic/particleFlux"]
val1 = val1[...]
val1_1d = np.transpose(val1[:, 0, 0, 0])
y1[0] = val1_1d[t]

val2 = f["/Diagnostic/particleFlux"]
val2 = val2[...]
val2_1d = np.transpose(val2[:, 0, 1, 0])
y2[0] = val2_1d[t]



##============ Case1 ==========

f=h5.File("IC1.0/data_global.h5")

val1 = f["/Diagnostic/particleFlux"]
val1 = val1[...]
val1_1d = np.transpose(val1[:, 0, 0, 0])
y1[1] = val1_1d[t]

val2 = f["/Diagnostic/particleFlux"]
val2 = val2[...]
val2_1d = np.transpose(val2[:, 0, 1, 0])
y2[1] = val2_1d[t]

##============ Case2 ===========
f=h5.File("IC2.0/data_global.h5")

val1 = f["/Diagnostic/particleFlux"]
val1 = val1[...]
val1_1d = np.transpose(val1[:, 0, 0, 0])
y1[2] = val1_1d[t]

val2 = f["/Diagnostic/particleFlux"]
val2 = val2[...]
val2_1d = np.transpose(val2[:, 0, 1, 0])
y2[2] = val2_1d[t]


##============ Case3 ===========
f=h5.File("IC3.0/data_global.h5")

val1 = f["/Diagnostic/particleFlux"]
val1 = val1[...]
val1_1d = np.transpose(val1[:, 0, 0, 0])
y1[3] = val1_1d[t]

val2 = f["/Diagnostic/particleFlux"]
val2 = val2[...]
val2_1d = np.transpose(val2[:, 0, 1, 0])
y2[3] = val2_1d[t]


#pFlux = y1 + y2
pFlux = y2
pFlux_norm  = pFlux[0]	# used to normalize pFlux
print "pFlux: ", pFlux_norm
pFlux = pFlux / pFlux_norm




sp_temp1.yaxis.set_major_formatter(yformatter)
cf_temp1=sp_temp1.plot(x, pFlux, marker = '8', label = r'$\mathrm{D^+}$ ion particle flux')
cf_temp2=sp_temp1.plot(x, hFlux, marker = '^', label = r'Total heat flux')

sp_temp1.set_ylim((0.75, 1.02))

major_ticks = np.arange(0.75, 1.02, 0.05)                                              
#minor_ticks = np.arange(0, 31, 5)                                          
sp_temp1.set_yticks(major_ticks)                                                       
#sp_temp1.set_yticks(minor_ticks, minor=True)  




sp_temp1.grid(True)
sp_temp1.legend(fontsize = legend_fontsize, loc = 1)

sp_temp1.xaxis.set_ticks(x)

#sp_temp1.set_xlim((xmin, xmax))
#sp_temp1.set_ylim((ymin, ymax))

sp_temp1.set_xlabel(r'$\Gamma_\mathrm{C}\ \mathrm{(m^{-2}s^{-1})}$', fontsize = label_fontsize)
sp_temp1.set_ylabel(r'Normalized fluxes', fontsize = label_fontsize)





fig.savefig("all_Flux_C.pdf", dpi = 300)
##fig.show()       #when the program finishes,the figure disappears
#plt.axis('equal')
#plt.show()         #The command is OK


##1d plot=============================
##ax1=fig.add_subplot(2,1,1)
#flux=f["/1d/pflux"]
#flux=flux[...]

#line1,=ax1.plot(flux[0,:],flux[1,:])
#line2,=ax1.plot(flux[0,:],flux[2,:])
