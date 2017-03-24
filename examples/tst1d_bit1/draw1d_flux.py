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
import math
import ConfigParser



import h5py as h5
import numpy as np
import matplotlib.pyplot as plt


##inite the fig of matplotlib
fig=plt.figure(figsize=(10,8))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.4)
t = 14


##read data from file
f=h5.File("data_global.h5")
print f.keys()

group = f['/Fields']


val1 = f["/Diagnostic/particleNumber"]
val1 = val1[...]

nx = (val1.shape)[0]
x = np.linspace(0,nx,nx)
print nx


#============Total particle number======================================= 

val1 = f["/Diagnostic/particleNumber"]
val1 = val1[...]
val1_1d = np.transpose(val1[:, 0, 0, 0])

val2 = f["/Diagnostic/particleNumber"]
val2 = val2[...]
val2_1d = np.transpose(val2[:, 0, 0, 1])

sp_temp1=fig.add_subplot(3,1,1)

cf_temp1=sp_temp1.plot(x, val1_1d, label='electron')
cf_temp1=sp_temp1.plot(x, val2_1d, label='D1 ion')

xmin = x.min()
xmax = 20 #x.max()
ymin = 1.2 * min( val1_1d.min(), val2_1d.min() )
ymax = 1.2 * max( val1_1d.max(), val2_1d.max() ) 


#sp_temp1.plot([xmin, xmax],[-Va_D, -Va_D])


sp_temp1.axis([xmin,xmax,ymin, ymax])
#sp_temp1.set_yticks(np.arange(0,y.max(),100))
sp_temp1.set_xlabel('time')
sp_temp1.set_ylabel('Particle number')
sp_temp1.legend()



##============Particle flux======================================================
val1 = f["/Diagnostic/particleFlux"]
val1 = val1[...]
val1_1d = np.transpose(val1[:, 0, 0, 0])

val2 = f["/Diagnostic/particleFlux"]
val2 = val2[...]
val2_1d = np.transpose(val2[:, 0, 1, 0])


sp_temp1=fig.add_subplot(3,1,2)


cf_temp1=sp_temp1.plot(x, val1_1d, label='electron')
cf_temp1=sp_temp1.plot(x, val2_1d, label='D1 ion')

ymin = 1.2 * min( val1_1d.min(), val2_1d.min() )
ymax = 1.2 * max( val1_1d.max(), val2_1d.max() ) 



sp_temp1.axis([xmin,xmax,ymin, ymax])
#sp_temp1.set_yticks(np.arange(0,y.max(),100))
sp_temp1.set_xlabel('time')
sp_temp1.set_ylabel('Particle Flux')
sp_temp1.legend()






##============Eenergy flux======================================================

val1 = f["/Diagnostic/heatFlux"]
val1 = val1[...]
val1_1d = np.transpose(val1[:, 0, 0, 0])

val2 = f["/Diagnostic/heatFlux"]
val2 = val2[...]
val2_1d = np.transpose(val2[:, 0, 1, 0])


sp_temp1=fig.add_subplot(3,1,3)


cf_temp1=sp_temp1.plot(x, val1_1d, label='electron')
cf_temp1=sp_temp1.plot(x, val2_1d, label='D1 ion')

ymin = 1.2 * min( val1_1d.min(), val2_1d.min() )
ymax = 1.2 * max( val1_1d.max(), val2_1d.max() ) 

sp_temp1.axis([xmin,xmax,ymin, ymax])
#sp_temp1.set_yticks(np.arange(0,y.max(),100))
sp_temp1.set_xlabel('time')
sp_temp1.set_ylabel('heatFlux')



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

