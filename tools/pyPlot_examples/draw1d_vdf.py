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
f=h5.File("data_global.h5")
print f.keys()


t0 = 1.0e-11
ntime_step_avg = 10000
dt = t0 * ntime_step_avg / 1.0e-6


##============ VDF ======================================================
sp_temp1 = fig.add_subplot(1,1,1)

val = f["/VDF/VDF_tot_C"]
val = val[...]

for t in range(1, 700 , 100):
	time = t * dt
	val_1d = np.transpose(val[t, 0, 0, :])
	nx = (val_1d.shape)[0]
	x = np.linspace(0.0, nx, nx)
	
	cf_temp1=sp_temp1.plot(x, val_1d, linewidth=2.0, label = str(time)+"$\mu$s" )
	

sp_temp1.set_xlabel('Velocity')
sp_temp1.set_title('C velocity density function')



plt.legend()
fig.savefig("fig.png", dpi = 300)
#plt.axis('equal')
plt.show()         #The command is OK



