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
fig=plt.figure(figsize=(10,8))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.4)
t = 1


##read data from file
f=h5.File("data_global.h5")
print f.keys()

group = f['/Fields']
dims = group.attrs['dims_global']
dims = dims[...]

print dims

nx = dims[3]


dx=1.0
x = np.linspace(0,100.0,nx)
print nx


##============rho======================================================
val = f["/Fields/Rho_global"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])


sp_temp1=fig.add_subplot(2,1,1)


cf_temp1=sp_temp1.plot(x, val_1d)



sp_temp1.axis([x.min(),x.max(),val_1d.min(),val_1d.max()])
#sp_temp1.set_yticks(np.arange(0,y.max(),100))
sp_temp1.set_xlabel('x(mm)')
sp_temp1.set_ylabel('Charge density')




##============potential======================================================

val = f["/Fields/Phi_global"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])


sp_temp1=fig.add_subplot(2,1,2)


cf_temp1=sp_temp1.plot(x, val_1d)



sp_temp1.axis([x.min(),x.max(),val_1d.min(),val_1d.max()])
#sp_temp1.set_yticks(np.arange(0,y.max(),100))
sp_temp1.set_xlabel('x(mm)')
sp_temp1.set_ylabel('Electric potential')



fig.savefig("fig.png")
##fig.show()       #when the program finishes,the figure disappears
#plt.axis('equal')
plt.show()         #The command is OK


##1d plot=============================
##ax1=fig.add_subplot(2,1,1)
#flux=f["/1d/pflux"]
#flux=flux[...]

#line1,=ax1.plot(flux[0,:],flux[1,:])
#line2,=ax1.plot(flux[0,:],flux[2,:])
