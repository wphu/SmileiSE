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
fig=plt.figure()
#fig=plt.figure(figsize=(10,8))
#fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.4)
t = 9


##read data from file
f=h5.File("Fields_global.h5")
print f.keys()

group = f['/Fields']
dims = group.attrs['dims_global']
dims = dims[...]

print dims

nx = dims[2]
ny = dims[3]

dx=1.0
dy=1.0

y,x=np.mgrid[slice(dy,dy*(ny+0.5),dy), slice(dx,dx*(nx+1),dx)]
#y,x=np.mgrid[slice(dx,dx*(nx+1),dx), slice(dy,dy*(ny+0.5),dy)]
print nx, ny






##============rho======================================================
val = f["/Fields/Rho_global"]
val = val[...]

val_2d = np.transpose(val[t, 0, :, :])
#val_2d = val[t, 0, :, :]
sp_temp1=fig.add_subplot(2,1,1) #, aspect='equal')
levels=MaxNLocator(nbins=100).tick_values(val_2d.min(),val_2d.max())

if(val_2d.min() == val_2d.max()):
	ticks_val=np.linspace(val_2d.min(),val_2d.max()+1.0,5)
else:
	ticks_val=np.linspace(val_2d.min(),val_2d.max(),5)

print ticks_val
cf_temp1=sp_temp1.contourf(x,y,val_2d,cmap=cm.get_cmap('jet'),levels=levels)
cbar_temp1=fig.colorbar(cf_temp1,ticks=ticks_val)



sp_temp1.axis([x.min(),x.max(),y.min(),y.max()])
sp_temp1.set_yticks(np.arange(0,y.max(),100))
sp_temp1.set_title('rho')
sp_temp1.set_xlabel('x(mm)')
sp_temp1.set_ylabel('y(mm)')




##============potential======================================================
val = f["/Fields/Phi_global"]
val = val[...]

val_2d = np.transpose(val[t, 0, :, :])
#val_2d = val[t, 0, :, :]
sp_temp1=fig.add_subplot(2,1,2) #, aspect='equal')
levels=MaxNLocator(nbins=100).tick_values(val_2d.min(),val_2d.max())

if(val_2d.min() == val_2d.max()):
	ticks_val=np.linspace(val_2d.min(),val_2d.max()+1.0,5)
else:
	ticks_val=np.linspace(val_2d.min(),val_2d.max(),5)

print val_2d.min(),val_2d.max()
print ticks_val
cf_temp1=sp_temp1.contourf(x,y,val_2d,cmap=cm.get_cmap('jet'),levels=levels)
cbar_temp1=fig.colorbar(cf_temp1,ticks=ticks_val)



sp_temp1.axis([x.min(),x.max(),y.min(),y.max()])
sp_temp1.set_yticks(np.arange(0,y.max(),100))
sp_temp1.set_title('Potential')
sp_temp1.set_xlabel('x(mm)')
sp_temp1.set_ylabel('y(mm)')





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

