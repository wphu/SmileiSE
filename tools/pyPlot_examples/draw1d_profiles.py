from template import *


t = 7

n0 = 1.0e19
Ex0 = 1.0e6

x_step = 20

amplification_factor = 80.0


##read data from file
f=h5.File("data_global.h5")

group = f['/Fields']
dims = group.attrs['dims_global']
dims = dims[...]
nx = dims[3]


dx = 0.5e-5  # unit (m)
x = np.linspace(0, nx * dx, nx)
x = x * amplification_factor
x_less = x[::x_step]


xmin = x.min()
xmax = x.max()

x0 = int(nx / 2)
x1 = 100




##inite the fig of matplotlib
fig=plt.figure(figsize=(10,8))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.55)

##============potential======================================================
ax0=fig.add_subplot(3,1,1)

val = f["/Fields/Phi_global_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
print( "potential max: ", val_1d.max() )
print( "potential: ",val_1d[x0], val_1d[x1] )
val_1d = val_1d[::x_step]

line0=ax0.plot(x_less, val_1d, label = r"$\phi$", color='#1f77b4')

ax0.set_ylabel(r"$\phi\ \mathrm{(V)}$", color='#1f77b4', fontsize = label_fontsize)
ax0.tick_params('y', colors='#1f77b4')

major_ticks = np.arange(0, 91, 30)
ax0.set_yticks(major_ticks)




#double y axis
ax0_twinx = ax0.twinx()
val = f["/Fields/Ex_global_avg"]
val = val[...] / Ex0
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]

ax0_twinx.yaxis.set_major_formatter(yformatter)
line0=ax0_twinx.plot(x_less, val_1d, label = r"$E_\mathrm{x}$", linestyle = linestyles[1] , color='#ff7f0e')

ax0_twinx.set_ylabel(r"$E_x \ \mathrm{(10^6V/m)}$", color='#ff7f0e', fontsize = label_fontsize)
ax0_twinx.tick_params('y', colors='#ff7f0e')

lines1, labels1 = ax0.get_legend_handles_labels()
lines2, labels2 = ax0_twinx.get_legend_handles_labels()
ax0_twinx.legend(lines1 + lines2, labels1 + labels2, loc = 1, framealpha=1)

ax0.grid(True)



ax0.set_xlim((xmin, xmax))
ax0.annotate(r"$\mathbf{(a)}$", xy=get_axis_limits(ax0), annotation_clip=False)

##============rho======================================================
ax0=fig.add_subplot(3,1,2)

val = f["/Fields/Rho_global_e_avg"]
val = val[...] / n0
val_1d = np.transpose(val[t, 0, 0, :])
print( "Electron density: ",val_1d[x0], val_1d[x1] )
val0_1d = val_1d
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r"Electron")




val = f["/Fields/Rho_global_D1_avg"]
val = val[...] / n0
val_1d = np.transpose(val[t, 0, 0, :])
val1_1d = val_1d
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$\mathrm{D^+}$ ion', linestyle = linestyles[1])

ax0.grid(True)

'''
x0_ax_inset0 = x.max() * 0.4
y0_ax_inset0 = val0_1d.max() * 0.27
width0_ax_inset0 = val0_1d.max() * 0.6
height0_ax_inset0 = val0_1d.max() * 0.7
ax_inset0 = inset_axes(ax0, width = "60%", height = "80%" , bbox_to_anchor = (x0_ax_inset0, y0_ax_inset0, width0_ax_inset0, height0_ax_inset0), bbox_transform = ax0.transData)

x_mm = x * 1.0e3 / amplification_factor
val0_1d = val0_1d
val1_1d = val1_1d
ax_inset0.plot(x_mm, val0_1d)
ax_inset0.plot(x_mm, val1_1d, linestyle = linestyles[1])

ax_inset0.set_xlim(0.0, x_mm.max() * 0.02)
ax_inset0.set_ylim(0.0, 0.15)
ax_inset0.set_xlabel(r"$x\ \mathrm{(mm)}$", fontsize = inset_label_fontsize)
ax_inset0.xaxis.labelpad = 0.5

for tick in ax_inset0.xaxis.get_major_ticks():
	tick.label.set_fontsize(inset_label_fontsize)
for tick in ax_inset0.yaxis.get_major_ticks():
	tick.label.set_fontsize(inset_label_fontsize)
'''



ax0.legend(loc = 1, framealpha=1)
ax0.set_xlim((xmin, xmax))
ax0.set_ylim((0.0, 1.1))
major_ticks = np.arange(0, 1.1, 0.5)
ax0.set_yticks(major_ticks)
ax0.set_ylabel(r"$n\ \mathrm{(10^{19}m^{-3})}$", fontsize = label_fontsize)

ax0.annotate(r"$\mathbf{(b)}$", xy=get_axis_limits(ax0), annotation_clip=False)

##============ Temperature ======================================================
ax0=fig.add_subplot(3,1,3)

val = f["/Fields/T_global_e_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
print( "Electron Temperature: ",val_1d[x0], val_1d[x1] )
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = "Electron")


ax0.grid(True)

val = f["/Fields/T_global_D1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
print( "D+1 ion temperature: ",val_1d[x0], val_1d[x1] )
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$\mathrm{D^+}$ ion', linestyle = linestyles[1])


ax0.legend(loc = 1, framealpha=1)


ymin = 0
ymax = val_1d.max() * 1.2
ax0.set_xlim((xmin, xmax))
ax0.set_ylim((ymin, 91))


major_ticks = np.arange(0, 91, 30)
ax0.set_yticks(major_ticks)

ax0.set_xlabel(r"$x\ \mathrm{(m)}$", fontsize = label_fontsize)
ax0.set_ylabel(r"$T\ \mathrm{(eV)}$", fontsize = label_fontsize)

ax0.annotate(r"$\mathbf{(c)}$", xy=get_axis_limits(ax0), annotation_clip=False)



fig.savefig("Profiles.pdf", dpi = 300)
##fig.show()       #when the program finishes,the figure disappears
#plt.axis('equal')
#plt.show()         #The command is OK
