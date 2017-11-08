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
    print ax.get_xlim()[0], ax.get_xlim()[1]
    return ax.get_xlim()[0] + (ax.get_xlim()[1] - ax.get_xlim()[0]) * x_scale, ax.get_ylim()[1] + ( ax.get_ylim()[1] - ax.get_ylim()[0] ) * (y_scale - 1.0)



##inite the fig of matplotlib
fig=plt.figure(figsize=(10,8))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.55)



##read data from file
f=h5.File("restore/Restore0_global.h5")
print f.keys()


v_number = 150
# for electron
m_ov_2T = 9.109382616e-31 / (2.0 * 64.5 * 1.602e-19)

# for D+ ion
#m_ov_2T = 2.0 * 1.67262158e-27 / (2.0 * 60 * 1.602e-19)

##============ velocity distribution in x direction ================================
sp_temp1 = fig.add_subplot(2,1,1)

val = f["/D1/momentum0"]
val = val[...]
particle_number = val.shape[0]
v = np.zeros(particle_number)



#======= vx from simulation
val = f["/D1/momentum0"]
v0 = val[...]
for iv in np.arange(0, particle_number):
	v[iv] = v0[iv]

vmin = v.min() * 1.0
vmax = -vmin
v_step = (vmax - vmin) / v_number
v_step_half = 0.5 * v_step

v_x = np.zeros(v_number)
v_y = np.zeros(v_number)
v_y_maxwell = np.zeros(v_number)

v_x[0] = vmin
for i in np.arange(1, v_number):
	v_x[i] = v_x[i-1] + v_step

for v_i in v:
	i = int( (v_i -vmin + v_step_half) / v_step )
	if i < 0:
		v_y[0] = v_y[0] + 1
	elif i >= v_number:
		v_y[v_number-1] = v_y[v_number-1] + 1
	else:
		v_y[i] = v_y[i] + 1
v_y = 1.0 * v_y / (particle_number * v_step)

cf_temp1=sp_temp1.plot(v_x, v_y, label = r'PIC-$v_x$')


#======= vy from simulation
val = f["/D1/momentum1"]
v0 = val[...]
for iv in np.arange(0, particle_number):
	v[iv] = v0[iv]



vmin = v.min() * 1.0
vmax = -vmin
v_step = (vmax - vmin) / v_number
v_step_half = 0.5 * v_step

v_x = np.zeros(v_number)
v_y = np.zeros(v_number)
v_y_maxwell = np.zeros(v_number)

v_x[0] = vmin
for i in np.arange(1, v_number):
	v_x[i] = v_x[i-1] + v_step

for v_i in v:
	i = int( (v_i -vmin + v_step_half) / v_step )
	if i < 0:
		v_y[0] = v_y[0] + 1
	elif i >= v_number:
		v_y[v_number-1] = v_y[v_number-1] + 1
	else:
		v_y[i] = v_y[i] + 1
v_y = 1.0 * v_y / (particle_number * v_step)

cf_temp1=sp_temp1.plot(v_x, v_y, label = r'PIC-$v_y$')


#======= vz from simulation
val = f["/D1/momentum2"]
v0 = val[...]
for iv in np.arange(0, particle_number):
	v[iv] = v0[iv]



vmin = v.min() * 1.0
vmax = -vmin
v_step = (vmax - vmin) / v_number
v_step_half = 0.5 * v_step

v_x = np.zeros(v_number)
v_y = np.zeros(v_number)
v_y_maxwell = np.zeros(v_number)

v_x[0] = vmin
for i in np.arange(1, v_number):
	v_x[i] = v_x[i-1] + v_step

for v_i in v:
	i = int( (v_i -vmin + v_step_half) / v_step )
	if i < 0:
		v_y[0] = v_y[0] + 1
	elif i >= v_number:
		v_y[v_number-1] = v_y[v_number-1] + 1
	else:
		v_y[i] = v_y[i] + 1
v_y = 1.0 * v_y / (particle_number * v_step)

cf_temp1=sp_temp1.plot(v_x, v_y, label = r'PIC-$v_z$')



#======= from theory
for i in np.arange(0, v_number):
	v_y_maxwell[i] = math.sqrt(m_ov_2T/3.14) * math.exp( -m_ov_2T * v_x[i] * v_x[i] )

cf_temp1=sp_temp1.plot(v_x, v_y_maxwell, label = r'Theory')

sp_temp1.set_xlim((vmin, vmax))
sp_temp1.set_ylim((0, 2.3e-7))

sp_temp1.grid()
sp_temp1.legend(framealpha=1)

sp_temp1.xaxis.set_major_formatter(yformatter)
sp_temp1.yaxis.set_major_formatter(yformatter)

major_ticks = np.arange(0, 2.3e-7, 1.0e-7)                                                                                      
sp_temp1.set_yticks(major_ticks)   

sp_temp1.set_xlabel(r'Velocity (m/s)')
#sp_temp1.set_ylabel(r'$F_M \ (v_x)$')
sp_temp1.set_ylabel(r'$F_M$')

sp_temp1.annotate(r"$\mathbf{(a)}$", xy=get_axis_limits(sp_temp1), annotation_clip=False)


v_number = 150
##==================================== Speed distribution =====================================
sp_temp1 = fig.add_subplot(2,1,2)

val = f["/D1/momentum0"]
val = val[...]
particle_number = val.shape[0]
v = np.zeros(particle_number)

val = f["/D1/momentum0"]
v0 = val[...]
for iv in np.arange(0, particle_number):
	v[iv] = v[iv] + v0[iv] * v0[iv]

val = f["/D1/momentum1"]
v0 = val[...]
for iv in np.arange(0, particle_number):
	v[iv] = v[iv] + v0[iv] * v0[iv]

val = f["/D1/momentum2"]
v0 = val[...]
for iv in np.arange(0, particle_number):
	v[iv] = v[iv] + v0[iv] * v0[iv]


for iv in np.arange(0, particle_number):
	v[iv] = math.sqrt( v[iv] )



vmin = v.min() * 1.0
vmax = v.max() * 1.0
v_step = (vmax - vmin) / v_number
v_step_half = 0.5 * v_step

v_x = np.zeros(v_number)
v_y = np.zeros(v_number)
v_y_maxwell = np.zeros(v_number)

v_x[0] = vmin
for i in np.arange(1, v_number):
	v_x[i] = v_x[i-1] + v_step

#======= from theory
for i in np.arange(0, v_number):
	v_y_maxwell[i] = math.pow(m_ov_2T/3.14, 1.5) * math.exp( -m_ov_2T * v_x[i] * v_x[i] ) * 4.0 * 3.14 * v_x[i] * v_x[i]

#======= from simulation
for v_i in v:
	i = int( (v_i -vmin + v_step_half) / v_step )
	if i < 0:
		v_y[0] = v_y[0] + 1
	elif i >= v_number:
		v_y[v_number-1] = v_y[v_number-1] + 1
	else:
		v_y[i] = v_y[i] + 1
v_y = 1.0 * v_y / (particle_number * v_step)

cf_temp1=sp_temp1.plot(v_x, v_y, label = r'PIC')
cf_temp1=sp_temp1.plot(v_x, v_y_maxwell, label = r'Theory')

sp_temp1.set_xlim((vmin, vmax))
sp_temp1.set_ylim((0, 3.2e-7))

sp_temp1.grid()
sp_temp1.legend(framealpha=1)

sp_temp1.xaxis.set_major_formatter(yformatter)
sp_temp1.yaxis.set_major_formatter(yformatter)

major_ticks = np.arange(0, 3.2e-7, 1.0e-7)                                                                                      
sp_temp1.set_yticks(major_ticks)   

sp_temp1.set_xlabel(r'Speed (m/s)')
sp_temp1.set_ylabel(r'$f_M$')

sp_temp1.annotate(r"$\mathbf{(b)}$", xy=get_axis_limits(sp_temp1), annotation_clip=False)


#plt.legend()
fig.savefig("vdf_vx_D1.pdf", dpi = 300)
#plt.axis('equal')
#plt.show()         #The command is OK
