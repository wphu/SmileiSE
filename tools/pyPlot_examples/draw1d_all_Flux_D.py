from template import *

label_fontsize = 23
legend_fontsize = 23


t = 19

##read data from file
f=h5.File("ref/data_global.h5")

group = f['/Fields']
dims = group.attrs['dims_global']
dims = dims[...]


nx = dims[3]


dx = 0.5e-2  # unit (mm)
x = np.zeros(5)
y1 = np.zeros(5)
y2 = np.zeros(5)
y3 = np.zeros(5)
pFlux = np.zeros(5)		# particle flux
hFlux = np.zeros(5)		# heat flux

x[0] = 0.0
x[1] = 0.2
x[2] = 0.4
x[3] = 0.6
x[4] = 0.8


#xmin = x.min()
#xmax = x.max()
ymin = 0.0
ymax = 6.4e6

##inite the fig of matplotlib
fig=plt.figure(figsize=(10,8))
fig.subplots_adjust(top=0.9,bottom=0.1,wspace=0.5,hspace=0.55)

sp_temp1=fig.add_subplot(1,1,1)
# ================================= heat flux =========================================
##============Ref =======
f=h5.File("ref/data_global.h5")

val1 = f["/Diagnostic/heatFlux"]
val1 = val1[...]
val1_1d = np.transpose(val1[:, 0, 0, 0])
y1[0] = val1_1d[t]

val2 = f["/Diagnostic/heatFlux"]
val2 = val2[...]
val2_1d = np.transpose(val2[:, 0, 1, 0])
y2[0] = val2_1d[t]



##============ Case1 =========

f=h5.File("Re0.2/data_global.h5")

val1 = f["/Diagnostic/heatFlux"]
val1 = val1[...]
val1_1d = np.transpose(val1[:, 0, 0, 0])
y1[1] = val1_1d[13]

val2 = f["/Diagnostic/heatFlux"]
val2 = val2[...]
val2_1d = np.transpose(val2[:, 0, 1, 0])
y2[1] = val2_1d[13]

##============ Case2 ============
f=h5.File("Re0.4/data_global.h5")

val1 = f["/Diagnostic/heatFlux"]
val1 = val1[...]
val1_1d = np.transpose(val1[:, 0, 0, 0])
y1[2] = val1_1d[13]

val2 = f["/Diagnostic/heatFlux"]
val2 = val2[...]
val2_1d = np.transpose(val2[:, 0, 1, 0])
y2[2] = val2_1d[13]


##============ Case3 ===========
f=h5.File("Re0.6/data_global.h5")

val1 = f["/Diagnostic/heatFlux"]
val1 = val1[...]
val1_1d = np.transpose(val1[:, 0, 0, 0])
y1[3] = val1_1d[t]

val2 = f["/Diagnostic/heatFlux"]
val2 = val2[...]
val2_1d = np.transpose(val2[:, 0, 1, 0])
y2[3] = val2_1d[t]


##============ Case4 ===========
f=h5.File("Re0.8/data_global.h5")

val1 = f["/Diagnostic/heatFlux"]
val1 = val1[...]
val1_1d = np.transpose(val1[:, 0, 0, 0])
y1[4] = val1_1d[t]

val2 = f["/Diagnostic/heatFlux"]
val2 = val2[...]
val2_1d = np.transpose(val2[:, 0, 1, 0])
y2[4] = val2_1d[t]



hFlux = y1 + y2
hFlux_norm  = hFlux[0]	# used to normalize hFlux
print( "hFlux: ", hFlux_norm )
hFlux = hFlux / hFlux_norm




# ========================================== particle flux =======================
##============Ref =============
f=h5.File("ref/data_global.h5")

val1 = f["/Diagnostic/particleFlux"]
val1 = val1[...]
val1_1d = np.transpose(val1[:, 0, 0, 0])
y1[0] = val1_1d[t]

val2 = f["/Diagnostic/particleFlux"]
val2 = val2[...]
val2_1d = np.transpose(val2[:, 0, 1, 0])
y2[0] = val2_1d[t]



##============ Case1 ==========

f=h5.File("Re0.2/data_global.h5")

val1 = f["/Diagnostic/particleFlux"]
val1 = val1[...]
val1_1d = np.transpose(val1[:, 0, 0, 0])
y1[1] = val1_1d[13]

val2 = f["/Diagnostic/particleFlux"]
val2 = val2[...]
val2_1d = np.transpose(val2[:, 0, 1, 0])
y2[1] = val2_1d[13]

##============ Case2 ===========
f=h5.File("Re0.4/data_global.h5")

val1 = f["/Diagnostic/particleFlux"]
val1 = val1[...]
val1_1d = np.transpose(val1[:, 0, 0, 0])
y1[2] = val1_1d[13]

val2 = f["/Diagnostic/particleFlux"]
val2 = val2[...]
val2_1d = np.transpose(val2[:, 0, 1, 0])
y2[2] = val2_1d[13]


##============ Case3 ===========
f=h5.File("Re0.6/data_global.h5")

val1 = f["/Diagnostic/particleFlux"]
val1 = val1[...]
val1_1d = np.transpose(val1[:, 0, 0, 0])
y1[3] = val1_1d[t]

val2 = f["/Diagnostic/particleFlux"]
val2 = val2[...]
val2_1d = np.transpose(val2[:, 0, 1, 0])
y2[3] = val2_1d[t]

##============ Case4 ===========
f=h5.File("Re0.8/data_global.h5")

val1 = f["/Diagnostic/particleFlux"]
val1 = val1[...]
val1_1d = np.transpose(val1[:, 0, 0, 0])
y1[4] = val1_1d[t]

val2 = f["/Diagnostic/particleFlux"]
val2 = val2[...]
val2_1d = np.transpose(val2[:, 0, 1, 0])
y2[4] = val2_1d[t]




#pFlux = y1 + y2
pFlux = y2
pFlux_norm  = pFlux[0]	# used to normalize pFlux
print( "pFlux: ", pFlux_norm )
pFlux = pFlux / pFlux_norm




sp_temp1.yaxis.set_major_formatter(yformatter)
cf_temp1=sp_temp1.plot(x, pFlux, marker = '8', label = r'$\mathrm{D^+}$ ion particle flux', linestyle = linestyles[0])
cf_temp2=sp_temp1.plot(x, hFlux, marker = '^', label = r'Total heat flux', linestyle = linestyles[1])

sp_temp1.set_ylim((0, 1.02))

major_ticks = np.arange(0, 1.02, 0.2)
#minor_ticks = np.arange(0, 31, 5)
sp_temp1.set_yticks(major_ticks)
#sp_temp1.set_yticks(minor_ticks, minor=True)




sp_temp1.grid(True)
sp_temp1.legend(fontsize = legend_fontsize, loc = 1)

sp_temp1.xaxis.set_ticks(x)

#sp_temp1.set_xlim((xmin, xmax))
#sp_temp1.set_ylim((ymin, ymax))

sp_temp1.set_xlabel(r'$c_r$', fontsize = label_fontsize)
sp_temp1.set_ylabel(r'Normalized fluxes', fontsize = label_fontsize)





fig.savefig("all_Flux_D.pdf", dpi = 300)
##fig.show()       #when the program finishes,the figure disappears
#plt.axis('equal')
#plt.show()         #The command is OK


##1d plot=============================
##ax1=fig.add_subplot(2,1,1)
#flux=f["/1d/pflux"]
#flux=flux[...]

#line1,=ax1.plot(flux[0,:],flux[1,:])
#line2,=ax1.plot(flux[0,:],flux[2,:])
