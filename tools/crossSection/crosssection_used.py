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

#mpl.rcParams['text.usetex'] = True
#mpl.rcParams['text.latex.unicode'] = True
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

def get_axis_limits(ax, x_scale=0, y_scale=1.09):
    return ax.get_xlim()[1]*x_scale, ax.get_ylim()[1]*y_scale


##inite the fig of matplotlib
fig=plt.figure(figsize=(12,15))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.4)
t = 15




##============ D ionization and excitation ===========================================
ax1=fig.add_subplot(3,2,1)

labels = [r'Ion: $D-D^+$',
		  r'Exc: D 1s-2p',
		  r'Exc: D 1s-3p',
		  r'Exc: D 1s-4p',
		  r'Exc: D 1s-5p',
		  r'Exc: D 1s-6p',
		  r'Exc: D 1s-7p',
		  r'Exc: D 1s-8p',
		  r'Exc: D 1s-9p',
		  r'Exc: D 1s-10p',
		  r'CX:  $D-D^+$',
		 ]

filenames = ['H&D/data/Ionization_D_to_D+1.dat',
			 'H&D/data/Excitation_H_1s-2p.dat',
			 'H&D/data/Excitation_H_1s-3p.dat',
			 'H&D/data/Excitation_H_1s-4p.dat',
			 'H&D/data/Excitation_H_1s-5p.dat',
			 'H&D/data/Excitation_H_1s-6p.dat',
			 'H&D/data/Excitation_H_1s-7p.dat',
			 'H&D/data/Excitation_H_1s-8p.dat',
			 'H&D/data/Excitation_H_1s-9p.dat',
			 'H&D/data/Excitation_H_1s-10p.dat',
			 'H&D/data/Charge_exchange_H+-H.dat',
			]

n = len(filenames)
for i in np.arange(0, n):
	data1 = np.loadtxt(filenames[i])
	if(i == n-1):
		data1[:,0] = data1[:,0] * 2.0
		ax1_double = ax1.twinx()
		ax1_double.plot(data1[:,0], data1[:,1], label = labels[i], color = 'blue')
		ax1_double.set_ylabel(r'CX cross section', color = 'blue')
	else:
		ax1.plot(data1[:,0], data1[:,1], label = labels[i])


ax1.set_xlim((0, 500))
#ax1.set_ylabel(r'$T_e \ (eV)$')
#ax1.set_ylabel("Cross Section $(m^2)$")

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax1_double.get_legend_handles_labels()
ax1_double.legend(lines1 + lines2, labels1 + labels2, fontsize = 8)
#ax1.legend(fontsize = 8)
ax1.grid(True)
ax1.annotate('(a) Ion, Exc and CX of D', xy=get_axis_limits(ax1), annotation_clip=False)



##============ C ionization ===========================================
ax1=fig.add_subplot(3,2,2)

labels = [r'Ion: $D-D^+$',
		  r'Ion: $D^+-D^{+2}$',
		  r'Ion: $D^{+2}-D^{+3}$',
		 ]

filenames = ['C/data/Ionization_C_to_C+1.dat',
			 'C/data/Ionization_C+1_to_C+2.dat',
			 'C/data/Ionization_C+2_to_C+3.dat',
			]

n = len(filenames)
for i in np.arange(0, n):
	data1 = np.loadtxt(filenames[i])
	ax1.plot(data1[:,0], data1[:,1], label = labels[i])

ax1.set_xlim((0, 500))
#ax1.set_ylabel(r'$T_e \ (eV)$')
#ax1.set_ylabel("Cross Section $(m^2)$")

ax1.legend(fontsize = 8)
ax1.grid(True)
ax1.annotate('(b) Ionization of C', xy=get_axis_limits(ax1), annotation_clip=False)



##============ C excitation ===========================================
ax1=fig.add_subplot(3,2,3)


filenames = ['C/data/Excitation_C_2s22p21D-2s2p33D.dat',
			 'C/data/Excitation_C_2s22p21D-2s22p3s1P.dat',
			 'C/data/Excitation_C_2s22p23P-2s2p33D.dat',
			 'C/data/Excitation_C_2s22p23P-2s2p33P.dat',
			 'C/data/Excitation_C_2s22p23P-2s2p33S.dat',
			 'C/data/Excitation_C_2s22p23P-2s2p35S.dat',
			 'C/data/Excitation_C_2s22p23P-2s22p3d3D.dat',
			 'C/data/Excitation_C_2s22p23P-2s22p3p3P.dat',
			 'C/data/Excitation_C_2s22p23P-2s22p3s3P.dat',
			 'C/data/Excitation_C_2s22p23P-2s22p21D.dat',
			]

n = len(filenames)
for i in np.arange(0, n):
	data1 = np.loadtxt(filenames[i])
	ax1.plot(data1[:,0], data1[:,1])

ax1.set_xlim((0, 200))
#ax1.set_ylabel(r'$T_e \ (eV)$')
ax1.set_ylabel("Cross Section $(m^2)$", fontsize = 22)

#ax1.legend(fontsize = 8)
ax1.grid(True)
ax1.annotate('(c) Excitation of C', xy=get_axis_limits(ax1), annotation_clip=False)



##============ C+1 excitation ===========================================
ax1=fig.add_subplot(3,2,4)

filenames = ['C/data/Excitation_C+1_2s22p2P-2s2p22D.dat',
			 'C/data/Excitation_C+1_2s22p2P-2s2p22P.dat',
			 'C/data/Excitation_C+1_2s22p2P-2s2p22S.dat',
			 'C/data/Excitation_C+1_2s22p2P-2s2p24P.dat',
			 'C/data/Excitation_C+1_2s22p2P-2s23d2D.dat',
			 'C/data/Excitation_C+1_2s22p2P-2s23p2P.dat',
			 'C/data/Excitation_C+1_2s22p2P-2s23s2S.dat',
			]

n = len(filenames)
for i in np.arange(0, n):
	data1 = np.loadtxt(filenames[i])
	ax1.plot(data1[:,0], data1[:,1])

ax1.set_xlim((0, 200))
#ax1.set_ylabel(r'$T_e \ (eV)$')
#ax1.set_ylabel("Cross Section $(m^2)$")

#ax1.legend(fontsize = 8)
ax1.grid(True)
ax1.annotate(r'(d) Excitation of $C^+$', xy=get_axis_limits(ax1), annotation_clip=False)



##============ C+2 excitation ===========================================
ax1=fig.add_subplot(3,2,5)


filenames = ['C/data/Excitation_C+2_2p23P-2p21D.dat',
			 'C/data/Excitation_C+2_2s2p1P-2p21D.dat',
			 'C/data/Excitation_C+2_2s2p3P-2p23P.dat',
			]

n = len(filenames)
for i in np.arange(0, n):
	data1 = np.loadtxt(filenames[i])
	ax1.plot(data1[:,0], data1[:,1])

ax1.set_xlim((0, 200))
#ax1.set_ylabel(r'$T_e \ (eV)$')
#ax1.set_ylabel("Cross Section $(m^2)$")

#ax1.legend(fontsize = 8)
ax1.grid(True)
ax1.annotate(r'(e) Excitation of $C^{+2}$', xy=get_axis_limits(ax1), annotation_clip=False)


##============ C+3 excitation ===========================================
ax1=fig.add_subplot(3,2,6)


filenames = ['C/data/Excitation_C+3_2p2P-3d2D.dat',
			 'C/data/Excitation_C+3_2p2P-3s2S.dat',
			 'C/data/Excitation_C+3_2s2S-2p2P.dat',
			]

n = len(filenames)
for i in np.arange(0, n):
	data1 = np.loadtxt(filenames[i])
	ax1.plot(data1[:,0], data1[:,1])

ax1.set_xlim((0, 200))
#ax1.set_ylabel(r'$T_e \ (eV)$')
#ax1.set_ylabel("Cross Section $(m^2)$")

#ax1.legend(fontsize = 8)
ax1.grid(True)
ax1.annotate('(f) Excitation of $C^{+3}$', xy=get_axis_limits(ax1), annotation_clip=False)

ax1.annotate('$T_e \ (eV)$', xy=get_axis_limits(ax1, -0.4, -0.3), fontsize = 22, annotation_clip=False)







fig.savefig("Profiles.pdf", dpi = 300)
plt.show()
