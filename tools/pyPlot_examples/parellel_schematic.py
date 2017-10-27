import matplotlib
matplotlib.use('Agg')


import matplotlib.pyplot as plt
from numpy.random import rand
import numpy as np
import matplotlib as mpl


from mpl_toolkits.axes_grid1 import make_axes_locatable

mpl.rcParams['font.size'] = 16
mpl.rcParams['axes.linewidth'] = 2.0
#mpl.rcParams['font.weight'] = 'bold'
#mpl.rcParams['pdf.fonttype'] = 42


mpl.rcParams['xtick.major.size'] = 2
mpl.rcParams['ytick.major.size'] = 2

mpl.rcParams['lines.linewidth'] = 2.0



#mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['font.family'] = 'sans-serif'
#mpl.rcParams['mathtext.default'] = 'regular'
#mpl.rcParams['pdf.fonttype'] = 3









def get_axis_limits(ax, x_scale=0, y_scale=1.03):
    return ax.get_xlim()[1]*x_scale, ax.get_ylim()[1]*y_scale

xmin = 0
xmax = 3
ymin = 0
ymax = 3

offset = 0.1

fig = plt.figure(1, figsize=(13, 6))

#=========================== One =============================================
ax1 = plt.subplot2grid((1,1), (0, 0), aspect = 'equal') #('141', aspect='equal')
divider = make_axes_locatable(ax1)

color = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8']
cores = ['core1', 'core2', 'core3', 'core4', 'core5', 'core6', 'core7', 'core8', 'core9']
'''
for color in ['red', 'green', 'blue']:
    n = 750
    x, y = 3 * rand(2, n)
    scale = 20.0 * rand(n)
    ax.scatter(x, y, c=color, s=scale, label=color,
               alpha=0.3, edgecolors='none')
'''

for xstart in np.arange(xmin, xmax):
    for ystart in np.arange(ymin, ymax):
        n = 50
	icolor = xstart * 3 + ystart
	x = xstart + 0.5 * offset + 0.9 * rand(1, n)
	y = ystart + 0.5 * offset + 0.9 * rand(1, n)
	scale = 15
	ax1.scatter(x, y, c=color[icolor], s=scale, label=cores[icolor], alpha=0.6)


ax1.axvline(x = 1, ymin = ymin, ymax = ymax, c="black",zorder=10, clip_on=True)
ax1.axvline(x = 2, ymin = ymin, ymax = ymax, c="black",zorder=10, clip_on=True)

ax1.axhline(y = 1, xmin = xmin, xmax = xmax, c="black",zorder=10, clip_on=True)
ax1.axhline(y = 2, xmin = xmin, xmax = xmax, c="black",zorder=10, clip_on=True)


ax1.set_xticks(np.arange(xmin, xmax, 0.2))
ax1.set_yticks(np.arange(ymin, ymax, 0.2))
ax1.grid()

ax1.set_xlim((xmin, xmax))
ax1.set_ylim((ymin, ymax))
ax1.legend(bbox_to_anchor=(1.01, 0.02, 0.102, 1.0), loc=2,
       ncol=1, borderaxespad=0., fontsize = 10)

plt.tick_params(
    left = 'off',
    bottom = 'off',
    labelleft = 'off',
    labelbottom='off')

ax1.annotate('(a)', xy=get_axis_limits(ax1), annotation_clip=False)

#=========================== Two =============================================
#ax2 = fig.add_subplot('142', aspect='equal')
ax2 = divider.append_axes("right", size = "100%", pad = 1.3)


offset = 0.1 * xmax
n = 50
icolor = 0
x = 0.5 * offset + 0.9 * xmax * rand(1, n)
y = 0.5 * offset + 0.9 * ymax * rand(1, n)
scale = 15
ax2.scatter(x, y, c=color[icolor], s=scale, label=cores[icolor], alpha=0.6)


ax2.set_xticks(np.arange(xmin, xmax, 0.2))
ax2.set_yticks(np.arange(ymin, ymax, 0.2))
ax2.grid()

ax2.set_xlim((xmin, xmax))
ax2.set_ylim((ymin, ymax))
ax2.legend(bbox_to_anchor=(0., 1.02, 1.0, 0.12), loc=4,
       ncol=1, borderaxespad=0., fontsize = 10)

plt.tick_params(
    left = 'off',
    bottom = 'off',
    labelleft = 'off',
    labelbottom='off')

ax2.annotate('+', xy=get_axis_limits(ax2, 1.02, 0.5), annotation_clip=False)
ax2.annotate('(b)', xy=get_axis_limits(ax2), annotation_clip=False)
#=========================== Three =============================================
ax3 = divider.append_axes("right", size = "100%", pad = 0.3)

offset = 0.1 * xmax
n = 50
icolor = 1
x = 0.5 * offset + 0.9 * xmax * rand(1, n)
y = 0.5 * offset + 0.9 * ymax * rand(1, n)
scale = 15
ax3.scatter(x, y, c=color[icolor], s=scale, label=cores[icolor], alpha=0.6)





ax3.set_xticks(np.arange(xmin, xmax, 0.2))
ax3.set_yticks(np.arange(ymin, ymax, 0.2))
ax3.grid()

ax3.set_xlim((xmin, xmax))
ax3.set_ylim((ymin, ymax))
ax3.legend(bbox_to_anchor=(0., 1.02, 1.0, 0.12), loc=4,
       ncol=1, borderaxespad=0., fontsize = 10)

plt.tick_params(
    left = 'off',
    bottom = 'off',
    labelleft = 'off',
    labelbottom='off')
ax3.annotate('...', xy=get_axis_limits(ax3, 1.02, 0.505), annotation_clip=False)

#=========================== Four =============================================
ax4 = divider.append_axes("right", size = "100%", pad = 0.3)

offset = 0.1 * xmax
n = 50
icolor = 8
x = 0.5 * offset + 0.9 * xmax * rand(1, n)
y = 0.5 * offset + 0.9 * ymax * rand(1, n)
scale = 15
ax4.scatter(x, y, c=color[icolor], s=scale, label=cores[icolor], alpha=0.6)





ax4.set_xticks(np.arange(xmin, xmax, 0.2))
ax4.set_yticks(np.arange(ymin, ymax, 0.2))
ax4.grid()

ax4.set_xlim((xmin, xmax))
ax4.set_ylim((ymin, ymax))
ax4.legend(bbox_to_anchor=(0., 1.02, 1.0, 0.12), loc=4,
       ncol=1, borderaxespad=0., fontsize = 10)

plt.tick_params(
    left = 'off',
    bottom = 'off',
    labelleft = 'off',
    labelbottom='off')


plt.savefig('parallel_schematic.pdf', dpi = 300)
#plt.show()
