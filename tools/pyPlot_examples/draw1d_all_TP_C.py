from template import *

legend_fontsize = 14

t = 19

x_step = 20

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
fig=plt.figure(figsize=(10,8))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.55)

# ===================================== Te ================================
ax0=fig.add_subplot(3,1,1)
##============Ref =============
f=h5.File("ref/data_global.h5")

val = f["/Fields/T_global_e_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = "Ref case", linestyle = linestyles[0])


##============ Case1 ==========
f=h5.File("IC1.0/data_global.h5")

val = f["/Fields/T_global_e_avg"]
val = val[...]
val_1d = np.transpose(val[13, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$\Gamma_C = \mathrm{1.0\times 10^{22}m^{-2}s^{-1}}$', linestyle = linestyles[1])


##============ Case2 ============
f=h5.File("IC2.0/data_global.h5")

val = f["/Fields/T_global_e_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$\Gamma_C = \mathrm{2.0\times 10^{22}m^{-2}s^{-1}}$', linestyle = linestyles[2])



##============ Case3 ===========
f=h5.File("IC3.0/data_global.h5")

val = f["/Fields/T_global_e_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$\Gamma_C = \mathrm{3.0\times 10^{22}m^{-2}s^{-1}}$', linestyle = linestyles[3])



ax0.grid(True)
ax0.legend(loc = 1, framealpha=1, fontsize = legend_fontsize)
ax0.set_xlim((xmin, xmax))
ax0.set_ylim((0.0, 30.0))

major_ticks = np.arange(0, 31, 10)
ax0.set_yticks(major_ticks)

ax0.set_ylabel(r"$T_\mathrm{e} \ \mathrm{(eV)}$", fontsize = label_fontsize)


ax0.annotate(r"$\mathbf{(a)}$", xy=get_axis_limits(ax0), annotation_clip=False)


# ===================================== Ti ================================
ax0=fig.add_subplot(3,1,2)
ax0.yaxis.set_major_formatter(yformatter)
##============Ref =============
f=h5.File("ref/data_global.h5")

val = f["/Fields/T_global_D1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = "Ref case", linestyle = linestyles[0])


##============ Case1 ==========
f=h5.File("IC1.0/data_global.h5")

val = f["/Fields/T_global_D1_avg"]
val = val[...]
val_1d = np.transpose(val[13, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$\Gamma = 1.0\times 10^{22}m^{-2}s^{-1}$', linestyle = linestyles[1])


##============ Case2 ============
f=h5.File("IC2.0/data_global.h5")

val = f["/Fields/T_global_D1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$\Gamma = 2.0\times 10^{22}m^{-2}s^{-1}$', linestyle = linestyles[2])



##============ Case3 ===========
f=h5.File("IC3.0/data_global.h5")

val = f["/Fields/T_global_D1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$\Gamma = 3.0\times 10^{22}m^{-2}s^{-1}$', linestyle = linestyles[3])



ax0.grid(True)
#ax0.legend(loc = 1)
ax0.set_xlim((xmin, xmax))
ax0.set_ylim((0.0, 90.0))
ax0.set_ylabel(r"$T_{\mathrm{D^+}} \ \mathrm{(eV)}$", fontsize = label_fontsize)

major_ticks = np.arange(0, 91, 30)
ax0.set_yticks(major_ticks)

ax0.annotate(r"$\mathbf{(b)}$", xy=get_axis_limits(ax0), annotation_clip=False)




# =============================== Potential ================================
ax0=fig.add_subplot(3,1,3)
##============Ref ===============
f=h5.File("ref/data_global.h5")

val = f["/Fields/Phi_global_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = "Ref case", linestyle = linestyles[0])


##============ Case1 ============

f=h5.File("IC1.0/data_global.h5")

val = f["/Fields/Phi_global_avg"]
val = val[...]
val_1d = np.transpose(val[13, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = "InjectC 0.1", linestyle = linestyles[1])


##============ Case2 ============
f=h5.File("IC2.0/data_global.h5")

val = f["/Fields/Phi_global_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = "InjectC 0.5", linestyle = linestyles[2])



##============ Case3 ===========
f=h5.File("IC3.0/data_global.h5")

val = f["/Fields/Phi_global_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = "InjectC 1.0", linestyle = linestyles[3])



ax0.grid(True)
ax0.set_xlim((xmin, xmax))
ax0.set_ylim((0.0, 90.0))

major_ticks = np.arange(0, 91, 30)
ax0.set_yticks(major_ticks)


ax0.set_xlabel(r"$x\ \mathrm{(m)}$", fontsize = label_fontsize)
ax0.set_ylabel(r"$\phi\ \mathrm{(V)}$", fontsize = label_fontsize)


ax0.annotate(r"$\mathbf{(c)}$", xy=get_axis_limits(ax0), annotation_clip=False)


fig.savefig("all_TP_C.pdf", dpi = 300)
##fig.show()       #when the program finishes,the figure disappears
#plt.axis('equal')
#plt.show()         #The command is OK
