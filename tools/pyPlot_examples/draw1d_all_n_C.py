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


dx = 0.5e-5  # unit (m)
x = np.linspace(0, nx * dx, nx)

amplification_factor = 80.0
x = x * amplification_factor

xmin = x.min()
xmax = x.max()

x1=5000
x2=4500
x3=12





# ===================================== Ne ================================
sp_temp1=fig.add_subplot(3,1,1)
sp_temp1.yaxis.set_major_formatter(yformatter)
##============Ref =============
f=h5.File("ref/data_global.h5")

val = f["/Fields/Rho_global_e_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = r"Ref")
print "potential max: ", val_1d.max()
print "potential: ",val_1d[x1], val_1d[x2], val_1d[x3]

##============ Case1 ==========
f=h5.File("IC1.0/data_global.h5")

val = f["/Fields/Rho_global_e_avg"]
val = val[...]
val_1d = np.transpose(val[13, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = r'$\Gamma_C = \mathrm{1.0\times 10^{22}m^{-2}s^{-1}}$')
print "potential max: ", val_1d.max()
print "potential: ",val_1d[x1], val_1d[x2], val_1d[x3]

##============ Case2 ============
f=h5.File("IC2.0/data_global.h5")

val = f["/Fields/Rho_global_e_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = r'$\Gamma_C = \mathrm{2.0\times 10^{22}m^{-2}s^{-1}}$')
print "potential max: ", val_1d.max()
print "potential: ",val_1d[x1], val_1d[x2], val_1d[x3]


##============ Case3 ===========
f=h5.File("IC3.0/data_global.h5")

val = f["/Fields/Rho_global_e_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = r'$\Gamma_C = \mathrm{3.0\times 10^{22}m^{-2}s^{-1}}$')
print "potential max: ", val_1d.max()
print "potential: ",val_1d[x1], val_1d[x2], val_1d[x3]


sp_temp1.grid(True)
sp_temp1.legend(loc = 1, framealpha=1, fontsize = legend_fontsize)
sp_temp1.set_xlim((xmin, xmax))
#sp_temp1.set_ylim((0.0, 30.0))

major_ticks = np.arange(0, 1.2e19, 0.5e19)                                              
#minor_ticks = np.arange(0, 31, 5)                                          
sp_temp1.set_yticks(major_ticks)                                                       
#sp_temp1.set_yticks(minor_ticks, minor=True)  

sp_temp1.set_ylabel(r"$n_\mathrm{e} \ \mathrm{(m^{-3})}$", fontsize = label_fontsize)


sp_temp1.annotate(r"$\mathbf{(a)}$", xy=get_axis_limits(sp_temp1), annotation_clip=False)


# ===================================== Ni ================================
sp_temp1=fig.add_subplot(3,1,2)
sp_temp1.yaxis.set_major_formatter(yformatter)
##============Ref =============
f=h5.File("ref/data_global.h5")

val = f["/Fields/Rho_global_D1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = "Ref")
print "potential max: ", val_1d.max()
print "potential: ",val_1d[x1], val_1d[x2], val_1d[x3]

##============ Case1 ==========
f=h5.File("IC1.0/data_global.h5")

val = f["/Fields/Rho_global_D1_avg"]
val = val[...]
val_1d = np.transpose(val[13, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = r'$\Gamma = 1.0\times 10^{22}m^{-2}s^{-1}$')
print "potential max: ", val_1d.max()
print "potential: ",val_1d[x1], val_1d[x2], val_1d[x3]

##============ Case2 ============
f=h5.File("IC2.0/data_global.h5")

val = f["/Fields/Rho_global_D1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = r'$\Gamma = 2.0\times 10^{22}m^{-2}s^{-1}$')
print "potential max: ", val_1d.max()
print "potential: ",val_1d[x1], val_1d[x2], val_1d[x3]


##============ Case3 ===========
f=h5.File("IC3.0/data_global.h5")

val = f["/Fields/Rho_global_D1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = r'$\Gamma = 3.0\times 10^{22}m^{-2}s^{-1}$')
print "potential max: ", val_1d.max()
print "potential: ",val_1d[x1], val_1d[x2], val_1d[x3]


sp_temp1.grid(True)
#sp_temp1.legend(loc = 1)
sp_temp1.set_xlim((xmin, xmax))
#sp_temp1.set_ylim((0.0, 90.0))
sp_temp1.set_ylabel(r"$n_{\mathrm{D^+}} \ \mathrm{(m^{-3})}$", fontsize = label_fontsize)

major_ticks = np.arange(0, 1.2e19, 0.5e19)                                              
#minor_ticks = np.arange(0, 31, 5)                                          
sp_temp1.set_yticks(major_ticks)                                                       
#sp_temp1.set_yticks(minor_ticks, minor=True)  

sp_temp1.annotate(r"$\mathbf{(b)}$", xy=get_axis_limits(sp_temp1), annotation_clip=False)




# =============================== density of C ================================
sp_temp1=fig.add_subplot(3,1,3)
##============Ref ===============
f=h5.File("ref/data_global.h5")

val = f["/Fields/Rho_global_C_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = "Ref")
print "potential max: ", val_1d.max()
print "potential: ",val_1d[x1], val_1d[x2], val_1d[x3]

##============ Case1 ============

f=h5.File("IC1.0/data_global.h5")

val = f["/Fields/Rho_global_C_avg"]
val = val[...]
val_1d = np.transpose(val[13, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = "InjectC 0.1")
print "potential max: ", val_1d.max()
print "potential: ",val_1d[x1], val_1d[x2], val_1d[x3]

##============ Case2 ============
f=h5.File("IC2.0/data_global.h5")

val = f["/Fields/Rho_global_C_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = "InjectC 0.5")
print "potential max: ", val_1d.max()
print "potential: ",val_1d[x1], val_1d[x2], val_1d[x3]


##============ Case3 ===========
f=h5.File("IC3.0/data_global.h5")

val = f["/Fields/Rho_global_C_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
cf_temp1=sp_temp1.plot(x, val_1d, label = "InjectC 1.0")
print "potential max: ", val_1d.max()
print "potential: ",val_1d[x1], val_1d[x2], val_1d[x3]


sp_temp1.grid(True)
#sp_temp1.legend(loc = 1)
sp_temp1.set_xlim((xmin, xmax))
#sp_temp1.set_ylim((0.0, 90.0))

major_ticks = np.arange(0, 2.1e19, 0.5e19)                                              
#minor_ticks = np.arange(0, 31, 5)                                          
sp_temp1.set_yticks(major_ticks)                                                       
#sp_temp1.set_yticks(minor_ticks, minor=True)  


#legend1=sp_temp1.legend(loc=(.6,.76),fontsize=16)
#sp_temp1.axis([x.min(),x.max(),val_1d.min(),val_1d.max()])
#sp_temp1.set_yticks(np.arange(0,y.max(),100))
sp_temp1.set_xlabel(r"$x$ (m)", fontsize = label_fontsize)
sp_temp1.set_ylabel(r"$n_\mathrm{C} \ \mathrm{(m^{-3})}$", fontsize = label_fontsize)


sp_temp1.annotate(r"$\mathbf{(c)}$", xy=get_axis_limits(sp_temp1), annotation_clip=False)


fig.savefig("all_n_C.pdf", dpi = 300)
##fig.show()       #when the program finishes,the figure disappears
#plt.axis('equal')
#plt.show()         #The command is OK


##1d plot=============================
##ax1=fig.add_subplot(2,1,1)
#flux=f["/1d/pflux"]
#flux=flux[...]

#line1,=ax1.plot(flux[0,:],flux[1,:])
#line2,=ax1.plot(flux[0,:],flux[2,:])
