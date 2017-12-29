from template import *

t = 19
V0 = 1.0e6
V1 = 1.0e5

##read data from file
f=h5.File("data_global.h5")

group = f['/Fields']
dims = group.attrs['dims_global']
dims = dims[...]

nx = dims[3]


dx = 0.5e-2  # unit (mm)
x = np.linspace(0, nx * dx, nx)

xmin = x.min()
xmax = x.max()

# ion sound speed
Va1 = math.sqrt( 90.3 * 1.602e-19 / (2.0 * 1.67262158e-27) )
Va2 = math.sqrt( 120.0 * 1.602e-19 / (2.0 * 1.67262158e-27) )

##inite the fig of matplotlib
fig=plt.figure(figsize=(10,8))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.4)

##============ Parallel Velocity ======================================================
ax0=fig.add_subplot(2,1,1)
ax0.yaxis.set_major_formatter(yformatter)

val = f["/Fields/Vparallel_global_D1_avg"]
val = val[...] / V0
val_1d = np.transpose(val[t, 0, 0, :])
val1_1d = val_1d
line0=ax0.plot(x, val1_1d, label = r"$\mathrm{D^+}$ ion", linestyle = linestyles[0])

val = f["/Fields/Vparallel_global_e_avg"]
val = val[...] / V0
val_1d = np.transpose(val[t, 0, 0, :])
val0_1d = val_1d
line0=ax0.plot(x, val0_1d, label = r"Electron", linestyle = linestyles[1])



ax0.legend()
xmin = 0.0
xmax = 0.25
ax0.set_xlim((xmin, xmax))
ax0.set_xlim((xmin, xmax))
#ax0.set_yticks(np.arange(0,y.max(),100))
ax0.set_xlabel(r"$x\ \mathrm{(mm)}$", fontsize = label_fontsize)
ax0.set_ylabel(r"$V\ \mathrm{(10^6m/s)}$", fontsize = label_fontsize)

ax0.annotate(r"$\mathbf{(a)}$", xy=get_axis_limits(ax0), annotation_clip=False)

##============ Ion Parallel Velocity ======================================================
ax0=fig.add_subplot(2,1,2)
ax0.yaxis.set_major_formatter(yformatter)

val = f["/Fields/Vparallel_global_D1_avg"]
val = val[...] / V1
val_1d = np.transpose(val[t, 0, 0, :])
val1_1d = val_1d
line0=ax0.plot(x, val1_1d, label = "D+1", linestyle = linestyles[0])

Va1 = Va1 / V1
ax0.axhline(y = -Va1, xmin = 0.0, xmax = 1.0, c="red",zorder=0, clip_on=False, linestyle = linestyles[1])

xmin = 0.0
xmax = 0.25
ymin = val1_1d.min() * 1.1
ymax = 0.0
ax0.set_xlim((xmin, xmax))
ax0.set_ylim((ymin, ymax))

ax0.set_xlabel(r"$x\ \mathrm{(mm)}$", fontsize = label_fontsize)
ax0.set_ylabel(r"$V_{\mathrm{D^+}}\ \mathrm{(10^5m/s)}$", fontsize = label_fontsize)

major_ticks = np.arange(0, -1.2, -0.5)
ax0.set_yticks(major_ticks)

ax0.annotate(r"$\mathbf{(b)}$", xy=get_axis_limits(ax0), annotation_clip=False)




fig.savefig("velocity.pdf", dpi = 300)
##fig.show()       #when the program finishes,the figure disappears
#plt.axis('equal')
plt.show()         #The command is OK
