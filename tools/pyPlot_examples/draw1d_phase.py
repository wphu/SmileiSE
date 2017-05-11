##>>>The code is used to read data from hdf5 file
##>>>and plot on the screen and output figure file using matplotlib-python

import Tkinter as tk
from Tkinter import *

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib import cm
matplotlib.use('TkAgg')

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure

import numpy as np
from numpy import arange, sin, pi

import ConfigParser



import h5py as h5
import numpy as np
import matplotlib.pyplot as plt


##inite the fig of matplotlib
fig=plt.figure(figsize=(10,4))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.4)
t = 1


##read data from file
f=h5.File("restore/Restore0_global.h5")
print f.keys()


t0 = 1.0e-11
ntime_step_avg = 10000
dt = t0 * ntime_step_avg / 1.0e-6


##============ Phase plot in x direction ======================================================



sp_temp1 = fig.add_subplot(1,1,1)

val = f["/D1/position0"]
x = val[...]
xmin = x.min()
xmax = x.max()

val = f["/D1/momentum0"]
vx = val[...]
vxmin = vx.min()
vxmax = vx.max()

	
cf_temp1=sp_temp1.scatter(x, vx, s=3)
	
sp_temp1.set_xlim((xmin, xmax))
sp_temp1.set_ylim((vxmin, vxmax))

sp_temp1.set_xlabel('x')
sp_temp1.set_ylabel('vx')
sp_temp1.set_title('Phase')



#plt.legend()
fig.savefig("fig.png", dpi = 300)
#plt.axis('equal')
plt.show()         #The command is OK



