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

mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['font.family'] = 'sans-serif'
#mpl.rcParams['mathtext.default'] = 'regular'

mpl.rcParams['font.size'] = 16
mpl.rcParams['axes.linewidth'] = 2.0
#mpl.rcParams['font.weight'] = 'bold'

mpl.rcParams['xtick.major.size'] = 2
mpl.rcParams['ytick.major.size'] = 2

mpl.rcParams['lines.linewidth'] = 2.0

#mpl.rcParams['grid.linestyle'] = ":"
mpl.rcParams['grid.linestyle'] = ":"
mpl.rcParams['grid.color'] = "black"

def get_axis_limits(ax, x_scale=0, y_scale=1.02):
    return ax.get_xlim()[1]*x_scale, ax.get_ylim()[1]*y_scale


##inite the fig of matplotlib
fig=plt.figure(figsize=(10,4))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.4)



##read data from file
f=h5.File("restore/Restore12_global.h5")
print f.keys()


v_number = 150
m_ov_2T = 9.109382616e-31 / (2.0 * 40 * 1.602e-19)

##============ Phase plot in x direction ======================================================
sp_temp1 = fig.add_subplot(1,1,1)

val = f["/e/momentum0"]
val = val[...]
particle_number = val.shape[0]
v = np.zeros(particle_number)

val = f["/e/momentum0"]
v0 = val[...]
for iv in np.arange(0, particle_number):
	v[iv] = v[iv] + v0[iv] * v0[iv]

val = f["/e/momentum1"]
v0 = val[...]
for iv in np.arange(0, particle_number):
	v[iv] = v[iv] + v0[iv] * v0[iv]

val = f["/e/momentum2"]
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

cf_temp1=sp_temp1.plot(v_x, v_y)
cf_temp1=sp_temp1.plot(v_x, v_y_maxwell)

sp_temp1.set_xlim((vmin, vmax))
#sp_temp1.set_ylim((ymin, ymax))

sp_temp1.set_xlabel('Speed $(m/s)$')
sp_temp1.set_ylabel('f')




#plt.legend()
fig.savefig("fig.png", dpi = 300)
#plt.axis('equal')
plt.show()         #The command is OK
