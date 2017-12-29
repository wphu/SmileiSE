from template import *

t = 19

n0 = 1.0e19

x_step = 40

##read data from file
f=h5.File("ref/data_global.h5")

group = f['/Fields']
dims = group.attrs['dims_global']
dims = dims[...]

nx = dims[3]


dx = 0.5e-5  # unit (m)
x = np.linspace(0, nx * dx, nx)

amplification_factor = 80.0
x = x * amplification_factor
x_less = x[::x_step]

xmin = x.min()
xmax = x.max()

##inite the fig of matplotlib
fig=plt.figure(figsize=(10,6))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.4)

# ===================================== Ne ================================
ax0=fig.add_subplot(2,1,1)
##============Ref =============
f=h5.File("ref/data_global.h5")

val = f["/Fields/Rho_global_e_avg"]
val = val[...] / n0
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = "Ref case", linestyle = linestyles[0])


##============ Case1 ==========
f=h5.File("Re0.2/data_global.h5")

val = f["/Fields/Rho_global_e_avg"]
val = val[...] / n0
val_1d = np.transpose(val[13, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$c_r$ = 0.2', linestyle = linestyles[1])


##============ Case2 ============
f=h5.File("Re0.4/data_global.h5")

val = f["/Fields/Rho_global_e_avg"]
val = val[...] / n0
val_1d = np.transpose(val[13, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$c_r$ = 0.4', linestyle = linestyles[2])


##============ Case3 ===========
f=h5.File("Re0.6/data_global.h5")

val = f["/Fields/Rho_global_e_avg"]
val = val[...] / n0
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$c_r$ = 0.6', linestyle = linestyles[3])


##============ Case4 ===========
f=h5.File("Re0.8/data_global.h5")

val = f["/Fields/Rho_global_e_avg"]
val = val[...] / n0
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$c_r$ = 0.8', marker = markers[0])


ax0.grid(True)
ax0.legend(loc = 1, framealpha=1, fontsize = legend_fontsize)
ax0.set_xlim((xmin, xmax))
#ax0.set_ylim((0.0, 30.0))

major_ticks = np.arange(0, 2.1, 0.5)
ax0.set_yticks(major_ticks)

ax0.set_ylabel(r"$n_\mathrm{e} \ \mathrm{(10^{19}m^{-3})}$", fontsize = label_fontsize)


ax0.annotate(r"$\mathbf{(a)}$", xy=get_axis_limits(ax0), annotation_clip=False)



# =============================== nD ================================
ax0=fig.add_subplot(2,1,2)
##============Ref ===============
f=h5.File("ref/data_global.h5")

val = f["/Fields/Rho_global_D_avg"]
val = val[...] / n0
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = "Ref case", linestyle = linestyles[0])

##============ Case1 ============

f=h5.File("Re0.2/data_global.h5")

val = f["/Fields/Rho_global_D_avg"]
val = val[...] / n0
val_1d = np.transpose(val[13, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$c_r=0.2$', linestyle = linestyles[1])

##============ Case2 ============
f=h5.File("Re0.4/data_global.h5")

val = f["/Fields/Rho_global_D_avg"]
val = val[...] / n0
val_1d = np.transpose(val[13, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$c_r=0.4$', linestyle = linestyles[2])


##============ Case3 ===========
f=h5.File("Re0.6/data_global.h5")

val = f["/Fields/Rho_global_D_avg"]
val = val[...] / n0
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$c_r=0.6$', linestyle = linestyles[3])


##============ Case4 ===========
f=h5.File("Re0.8/data_global.h5")

val = f["/Fields/Rho_global_D_avg"]
val = val[...] / n0
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$c_r=0.8$', marker = markers[0])


ax0.grid(True)
#ax0.legend(loc = 1)
ax0.set_xlim((xmin, xmax))
#ax0.set_ylim((0.0, 90.0))

major_ticks = np.arange(0, 9.1, 3.0)
ax0.set_yticks(major_ticks)


ax0.set_xlabel(r"$x\ \mathrm{(m)}$", fontsize = label_fontsize)
ax0.set_ylabel(r"$n_\mathrm{D} \ \mathrm{(10^{19}m^{-3})}$", fontsize = label_fontsize)


ax0.annotate(r"$\mathbf{(b)}$", xy=get_axis_limits(ax0), annotation_clip=False)


fig.savefig("all_n_D.pdf", dpi = 300)
##fig.show()       #when the program finishes,the figure disappears
#plt.axis('equal')
#plt.show()         #The command is OK
