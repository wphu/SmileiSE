from template import *


t = 19

n0 = 1.0e19
Ex0 = 1.0e6

##read data from file
f=h5.File("data_global.h5")

group = f['/Fields']
dims = group.attrs['dims_global']
dims = dims[...]

nx = dims[3]


dx = 0.5e-2  # unit (mm)
x = np.linspace(0, nx * dx, nx)

xmin = x.min()
xmax = x.max() * 0.01
x_sheath = 0.1

##inite the fig of matplotlib
fig=plt.figure(figsize=(10,8))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.55)

##============potential======================================================
ax0=fig.add_subplot(3,1,1)

val = f["/Fields/Phi_global_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])

line0=ax0.plot(x, val_1d, label = r"$\phi$", color='#1f77b4', linestyle = linestyles[0])

ax0.set_ylabel(r"$\phi\ \mathrm{(V)}$", color='#1f77b4', fontsize = label_fontsize)
ax0.tick_params('y', colors='#1f77b4')

major_ticks = np.arange(0, 91, 30)
ax0.set_yticks(major_ticks)



#double y axis: Ex
ax0_twinx = ax0.twinx()
val = f["/Fields/Ex_global_avg"]
val = val[...] / Ex0
val_1d = np.transpose(val[t, 0, 0, :])

ax0_twinx.yaxis.set_major_formatter(yformatter)
line0=ax0_twinx.plot(x, val_1d, label = r"$E_\mathrm{x}$", color='#ff7f0e', linestyle = linestyles[1])

ax0_twinx.set_ylabel(r"$E_x \ \mathrm{(10^6V/m)}$", color='#ff7f0e', fontsize = label_fontsize)
ax0_twinx.tick_params('y', colors='#ff7f0e')


ax0.grid(True)

lines1, labels1 = ax0.get_legend_handles_labels()
lines2, labels2 = ax0_twinx.get_legend_handles_labels()
ax0_twinx.legend(lines1 + lines2, labels1 + labels2, loc = 1, framealpha=1)


ax0.set_xlim((xmin, xmax))



ax0.annotate(r"$\mathbf{(a)}$", xy=get_axis_limits(ax0), annotation_clip=False)

##============rho======================================================
ax0=fig.add_subplot(3,1,2)

val = f["/Fields/Rho_global_e_avg"]
val = val[...] / n0
val_1d = np.transpose(val[t, 0, 0, :])

ax0.yaxis.set_major_formatter(yformatter)
line0=ax0.plot(x, val_1d, label = "Electron", linestyle = linestyles[0])

ymin = 0.0
ymax = val_1d.max()

val = f["/Fields/Rho_global_D1_avg"]
val = val[...] / n0
val_1d = np.transpose(val[t, 0, 0, :])
line0=ax0.plot(x, val_1d, label = r'$\mathrm{D^+}$ ion', linestyle = linestyles[1])

ax0.grid(True)
ax0.legend(loc = 1, framealpha=1)
ax0.set_xlim((xmin, xmax))
ax0.set_ylim((ymin, ymax))

ax0.set_ylabel(r"$n\ \mathrm{(10^{19}m^{-3})}$", fontsize = label_fontsize)

major_ticks = np.arange(0, 1.01, 0.5)
ax0.set_yticks(major_ticks)


ax0.annotate(r"$\mathbf{(b)}$", xy=get_axis_limits(ax0), annotation_clip=False)

##============ Temperature ======================================================
ax0=fig.add_subplot(3,1,3)

val = f["/Fields/T_global_e_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
line0=ax0.plot(x, val_1d, label = "Electron", linestyle = linestyles[0])

val = f["/Fields/T_global_D1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
line0=ax0.plot(x, val_1d, label = r'$\mathrm{D^+}$ ion', linestyle = linestyles[1])

ymin = 0.0
ymax = val_1d.max()
# Plot a line
ax0.axvline(x = x_sheath, ymin = 0.0, ymax = 4.1, c="red",zorder=0, clip_on=False)


ax0.grid(True)
ax0.legend(loc = 1, framealpha=1)
ax0.set_xlim((xmin, xmax))
ax0.set_ylim((ymin, 90))

major_ticks = np.arange(0, 91, 30)
ax0.set_yticks(major_ticks)



ax0.set_xlabel(r"$x\ \mathrm{(mm)}$", fontsize = label_fontsize)
ax0.set_ylabel(r"$T\ \mathrm{(eV)}$", fontsize = label_fontsize)

ax0.annotate(r"$\mathbf{(c)}$", xy=get_axis_limits(ax0), annotation_clip=False)




fig.savefig("Sheath.pdf", dpi = 300)
##fig.show()       #when the program finishes,the figure disappears
#plt.axis('equal')
plt.show()         #The command is OK
