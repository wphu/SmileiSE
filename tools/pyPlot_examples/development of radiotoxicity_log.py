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


mpl.rcParams['font.size'] = 18
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



##inite the fig of matplotlib
fig=plt.figure(figsize=(10,8))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.55)

x = np.zeros(6)
y0 = np.zeros(6)
y1 = np.zeros(6)
y2 = np.zeros(6)

x[0] = 0.0
x[1] = 10.0
x[2] = 30.0
x[3] = 100.0
x[4] = 300.0
x[5] = 500.0

y0[0] = 4.5e11
y0[1] = 2.2e11
y0[2] = 1.9e11
y0[3] = 1.3e11
y0[4] = 9.0e10
y0[5] = 7.5e10

y1[0] = 9.0e10
y1[1] = 3.0e9
y1[2] = 2.8e8
y1[3] = 3.0e7
y1[4] = 9.5e6
y1[5] = 5.1e6

y2[0] = 5.45e6
y2[1] = 5.4e6
y2[2] = 5.35e6
y2[3] = 5.3e6
y2[4] = 5.25e6
y2[5] = 5.2e6



ax0 = fig.add_subplot(1,1,1)

line0 = ax0.semilogy(x, y0, label = 'Fission PWR', marker = 's')
line0 = ax0.semilogy(x, y1, label = 'Fusion SEAFP (model 2)', marker = '^')
line0 = ax0.semilogy(x, y2, label = 'Coal ash', marker = '*')

ax0.grid()
ax0.legend()

ax0.set_xlim((0, 500))
ax0.set_ylim((1.0e6, 1.0e12))

ax0.set_xlabel('Year', fontweight='bold')
ax0.set_ylabel('Ingestion Doses (Sv)', fontweight='bold')


fig.savefig("Development of radiotoxicity_log.pdf", dpi = 300)
##fig.show()       #when the program finishes,the figure disappears
