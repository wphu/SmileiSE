from template import *

t = 19

n0 = 1.0e19

x_step = 20

##read data from file
f=h5.File("data_global.h5")
print f.keys()

group = f['/Fields']
dims = group.attrs['dims_global']
dims = dims[...]

print dims

nx = dims[3]


dx = 0.5e-5  # unit (m)
x = np.linspace(0, nx * dx, nx)

amplification_factor = 80.0
x = x * amplification_factor
x_less = x[::x_step]

xmin = x.min()
xmax = x.max() * 0.2
ymin = 0.0

##inite the fig of matplotlib
fig=plt.figure(figsize=(10,8))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.5)

##============rho======================================================
ax0=fig.add_subplot(2,1,1)

val = f["/Fields/Rho_global_e_avg"]
val = val[...] / n0
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'Electron', linestyle = linestyles[0])

val = f["/Fields/Rho_global_D1_avg"]
val = val[...] / n0
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$\mathrm{D^+ \ ion}$', linestyle = linestyles[1])


val = f["/Fields/Rho_global_C_avg"]
val = val[...] / n0
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$\mathrm{C \ atom}$', linestyle = linestyles[2])



val = f["/Fields/Rho_global_C1_avg"]
val = val[...] / n0
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$\mathrm{C^+ \ ion}$', linestyle = linestyles[3])

val = f["/Fields/Rho_global_C2_avg"]
val = val[...] / n0
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$\mathrm{C^{+2} \ ion}$', marker = markers[0])

val = f["/Fields/Rho_global_C3_avg"]
val = val[...] / n0
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$\mathrm{C^{+3} \ ion}$', marker = markers[1])


ax0.grid(True)
ax0.legend(loc = 1, framealpha = 1.0, fontsize = legend_fontsize)
ax0.set_xlim((xmin, xmax))
ax0.set_ylim((ymin))

major_ticks = np.arange(0, 2.1, 0.5)                                              
#minor_ticks = np.arange(0, 31, 5)                                          
ax0.set_yticks(major_ticks)                                                       
#ax0.set_yticks(minor_ticks, minor=True)  

#ax0.set_yticks(np.arange(0,y.max(),100))
#ax0.set_xlabel('x(mm)')
ax0.set_ylabel(r"$n\ \mathrm{(10^{19}m^{-3})}$", fontsize = label_fontsize)

ax0.annotate(r"$\mathbf{(a)}$", xy=get_axis_limits(ax0), annotation_clip=False)

##============ Temperature ======================================================
ax0=fig.add_subplot(2,1,2)

val = f["/Fields/T_global_e_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'Electron', linestyle = linestyles[0])

val = f["/Fields/T_global_D1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$D^+ \ ion$', linestyle = linestyles[1])
ymax = val_1d.max() * 1.1


val = f["/Fields/T_global_C_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'C atom', linestyle = linestyles[2])


val = f["/Fields/T_global_C1_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$C^+ \ ion$', linestyle = linestyles[3])


val = f["/Fields/T_global_C2_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$C^{+2} \ ion$', marker = markers[0])


val = f["/Fields/T_global_C3_avg"]
val = val[...]
val_1d = np.transpose(val[t, 0, 0, :])
val_1d = val_1d[::x_step]
line0=ax0.plot(x_less, val_1d, label = r'$C^{+3} \ ion$', marker = markers[1])

ymin = 0

ax0.grid(True)
#ax0.legend(loc = 1)
ax0.set_xlim((xmin, xmax))
ax0.set_ylim((ymin, 66))

major_ticks = np.arange(0, 71.0, 20.0)                                              
#minor_ticks = np.arange(0, 31, 5)                                          
ax0.set_yticks(major_ticks)                                                       
#ax0.set_yticks(minor_ticks, minor=True)  


#ax0.set_yticks(np.arange(0,y.max(),100))
ax0.set_xlabel(r"$x\ \mathrm{(m)}$", fontsize = label_fontsize)
ax0.set_ylabel(r"$T\ \mathrm{(eV)}$", fontsize = label_fontsize)

ax0.annotate(r"$\mathbf{(b)}$", xy=get_axis_limits(ax0), annotation_clip=False)

fig.savefig("nT_C.pdf", dpi = 300)
##fig.show()       #when the program finishes,the figure disappears
#plt.axis('equal')
#plt.show()         #The command is OK


##1d plot=============================
##ax1=fig.add_subplot(2,1,1)
#flux=f["/1d/pflux"]
#flux=flux[...]

#line1,=ax1.plot(flux[0,:],flux[1,:])
#line2,=ax1.plot(flux[0,:],flux[2,:])
